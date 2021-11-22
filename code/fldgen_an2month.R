library(gtools)
library(devtools)
library(dplyr)
library(doParallel)
library(parallel)
library(foreach)
library(fldgen)


#' Harmonize coordinate locations of the alpha fractions and climate field data
#'
#' Match alpha fraction grid cell geographic locations with the coordinates from
#' the climate field data.  This will left join the alpha fractions with the
#' coordinate locations from the fldgen emulator data.
#'
#' @param frac Object read in from the alpha fractions RDS file
#' @param frac_coordinates Coordinates object found in the alpha fractions RDS
#' @param fld_coordinates Coordinates object as generate by fldgen outputs
#' @param var Variable name of the input (e.g., prAdjust)
#'
#' @return reindex alpha fractions to match land grid cell index locations
reindex_grid <- function(frac, frac_coordinates, fld_coordinates, var){

  # Silence package checks
  '%>%' <- 'column_index' <- NULL

  # Check inputs
  stopifnot(is.data.frame(frac_coordinates))
  stopifnot(is.data.frame(fld_coordinates))

  frac_dim <- dim(frac_coordinates)
  fld_dim  <- dim(fld_coordinates)

  # If the dimensions of the field file is greater than the fraction file it means that the
  # fraction file resolution is less than the field
  if(any(frac_dim < fld_dim)) stop('Resolution of the monthly fraction coordinates is insufficient')
  if(any(!fld_coordinates$lat %in% frac_coordinates$lat)) stop('fraction file is missing lat values')
  if(any(!fld_coordinates$lon %in% frac_coordinates$lon)) stop('fraction file is missing lon values')

  # If the fraction file has larger dimensions then it could be caused by the NAs, use the
  # fld coordinate information to subset the fraction file.
  if(any(dim(frac_coordinates) > dim(fld_coordinates))){

    fld_coordinates %>%
      dplyr::rename(flds_column_index = column_index) %>%
      dplyr::left_join(frac_coordinates %>%
                         dplyr::rename(frac_index = column_index),
                       by = c('lat', 'lon')) ->
      lat_lon_mapping


    # Modify the fraction object so that it contains the correct grid cells.
    frac <- frac[[var]][ , lat_lon_mapping$frac_index]

  } else {

    frac <- frac[[var]]

  }

  frac
}


generate_monthly_data <- function(fld_data, alpha_fractions, fld_coordinates, startyr) {

    # number of years in the data
    nyear <- nrow(fld_data)

    # annual values are averages, but we need totals, so multiply these values by 12
    fld_data <- fld_data * 12

    # number of grid cells
    ngrid <- ncol(fld_data)

    # random seed
    rngstate <- .Random.seed

    # creates a matrix where every column is a grid cell where rows are months (nyear * 12 in length)
    #  returned values are the monthly fraction of each annual total per grid cell
    monthly_data <- foreach::foreach(igrid = 1:ngrid, .combine='cbind') %dopar% {

        ## Use the same rng state for each grid cell.  This prevents us from having
        ## excessive spatial variation within a single month.
        .Random.seed <<- rngstate

        # generate a matrix of random vectors for the target grid cell from the
        #  Dirichlet distribution having a column for each month and row for each year
        monthly_fractions <- gtools::rdirichlet(nyear, alpha_fractions[,igrid]) # matrix[nyear, 12]

        # vector of annual fldgen values for each year for the target grid cell
        annual_totals <- fld_data[, igrid] # vector[nyear]

        # multiply matrix of [yr, month] by vector [yr] to disaggregate annual
        #  totals to monthly using the Dirichlet fractions
        monthly_totals <- monthly_fractions * annual_totals # matrix[nyear, 12]

        # first:  transpose col=month, row=year to col=year, row=month
        #   next: squash to a vector where values are (year_1:month_1, year_1:month_2,...,year_2:month_1,...)
        as.vector(t(monthly_totals))

    }

    ## We need to add the time codes as row names.  We now have 12x as many rows
    ## as we used to, so we need to do some replicating
    fld_time <- seq(startyr, startyr + nyear - 1)

    year_code <- rep(fld_time, rep(12, length(fld_time))) # Each time code gets

    # repeated 12 times in succession
    month_code <- seq(from=0, length.out=length(year_code)) %% 12 + 1
    row.names(monthly_data) <- sprintf('%04d%02d', year_code, month_code)

    # Return the monthly data and the coordinates as a list
    res <- list(data = monthly_data, coordinates = fld_coordinates)

    res

}


generate_yearly_climate_fields <- function(emulator_file, tgav_file, seed_value, ngrids,
                                           nyear, startyr, scenario, variable='Tgav') {

    # load emulator from file
    emu <- readRDS(emulator_file)

    # set random seed value for field reproducibility
    set.seed(seed_value)

    # generate residuals
    resids <- fldgen::generate.TP.resids(emu, ngrids)

    tgavdf <- read.csv(tgav_file) %>%
        filter(variable == variable)  %>%
        filter(scenario == scenario) %>%
        filter(year >= startyr) %>%
        filter(year < startyr + nyear) %>%
        arrange(year, .by_group = FALSE)

    fullgrids <- fldgen::generate.TP.fullgrids(emu, resids, tgavdf$value)

    fld_coordinates <- emu$griddataT$coord %>% as.data.frame()

    fld_coordinates$column_index <- seq(nrow(fld_coordinates))

    return_list <- list(fullgrids=fullgrids, fld_coordinates=fld_coordinates)

    return_list

}


downscale_yearly_climate_fields <- function(alpha_fractions_file, yearly_climate_data, startyr, pr_var='pr',
                                            tas_var='tas') {

    # activate cluster for foreach library
    n_cores <- detectCores()
    doParallel::registerDoParallel(n_cores)

    print(paste0("Using ", n_cores, " cores."))

    print("Extracting yearly grids...")
    target_field <- yearly_climate_data$fullgrids$fullgrids[[1]]

    # read in alpha fraction file produced in an2month calibration
    print("Reading alpha fractions file...")
    alpha <- readRDS(alpha_fractions_file)

    # downscale precipitation
    print("Reindexing alpha grid for pr...")
    frac <- reindex_grid(frac=alpha,
                         frac_coordinates=alpha$coord,
                         fld_coordinates=yearly_climate_data$fld_coordinates,
                         var=pr_var)

    print("Downscaling variable:  pr...")
    t0 <- Sys.time()
    pr_result <- generate_monthly_data(fld_data=target_field$pr,
                                       alpha_fractions=frac,
                                       fld_coordinates=yearly_climate_data$fld_coordinates,
                                       startyr=startyr)
    print(paste0("Downscaled 'pr' in ", (Sys.time() - t0), " minutes."))

    # downscale tas
    print("Reindexing alpha grid for tas")
    frac <- reindex_grid(frac=alpha,
                         frac_coordinates=alpha$coord,
                         fld_coordinates=yearly_climate_data$fld_coordinates,
                         var=tas_var)


    print("Downscaling variable:  tas...")
    t0 <- Sys.time()
    tas_result <- generate_monthly_data(fld_data=target_field$tas,
                                        alpha_fractions=frac,
                                        fld_coordinates=yearly_climate_data$fld_coordinates,
                                        startyr=startyr)
    print(paste0("Downscaled 'tas' in ", (Sys.time() - t0), " minutes."))

    return(list('tas'=tas_result, 'pr'=pr_result))

}


generate_monthly_climate_fields <- function(emulator_file,
                                    tgav_file,
                                    seed_value,
                                    ngrids,
                                    nyear,
                                    startyr,
                                    scenario,
                                    alpha_fractions_file,
                                    variable='Tgav') {

    # generate yearly climate fields using fldgen
    print("Generating yearly climate fields...")
    t0 <- Sys.time()
    yearly_climate_data <- generate_yearly_climate_fields(emulator_file,
                            tgav_file, seed_value, ngrids, nyear, startyr,
                            scenario, variable='Tgav')
    print(paste0("Generated yearly climate fields in ", (Sys.time() - t0), " minutes."))

    # use an2month algorithm to downscale yearly climate fields to monthly
    print("Downscaling yearly climate fields to monthly resolution...")
    t0 <- Sys.time()
    results <- downscale_yearly_climate_fields(alpha_fractions_file,
                yearly_climate_data, startyr)
    print(paste0("Completed downscaling in ", (Sys.time() - t0), " minutes."))

    return(results)
}

