import gc
import os
import sys
import calendar

import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from xanthos import Xanthos

robjects.numpy2ri.activate()


def get_realizations(n_realizations, batch_size, slurm_idx):
    """Generate a list of realizations to process for the current node."""

    batches = range(0, n_realizations, batch_size)

    try:
        return [list(range(i, i + batch_size, 1)) for i in batches][slurm_idx]
    except IndexError:
        msg = f'No realization batch exists for SLURM_ARRAY_TASK_ID = {slurm_idx}.'
        raise IndexError(msg)


def collect_garbage():
    """Trigger manual garbage collection from R and Python."""

    r_garbage = robjects.r('gc()')
    py_garbage = gc.collect()


def generate_random_integer_seed():
    """Generate a random integer seed bounded by the min and max depth of np.int32.

    """

    # get a dytpe info object
    dtype_bounds = np.iinfo(np.int32)

    return np.random.randint(low=dtype_bounds.min, high=dtype_bounds.max)


def build_days_of_month_list(start_yr, n_years):
    """Generate a list of days of the month over a number of years that accounts for leap years."""

    days_per_month_leap = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    days_per_month_nonleap = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    days_per_month = []

    for i in range(start_yr, start_yr + n_years, 1):

        # add days of month to list based on leap status
        if calendar.isleap(i):
            days_per_month.extend(days_per_month_leap)
        else:
            days_per_month.extend(days_per_month_nonleap)

    return days_per_month


def convert_pr(arr, start_yr, n_years):
    """Convert precipitation from kg/m2/s to mm/month and transpose"""

    sec_in_min = 60
    min_in_hr = 60
    hr_in_day = 24

    # calculate the number of days per month for the year span accounting for leap years
    days_per_month = build_days_of_month_list(start_yr, n_years)

    return arr.T * sec_in_min * min_in_hr * hr_in_day * days_per_month


def generate_monthly_climate_fields(emulator_file,
                                    tgav_file,
                                    seed_value,
                                    ngrids,
                                    n_years,
                                    start_yr,
                                    scenario,
                                    alpha_fractions_file,
                                    r_script):
    """Runs fldgen and an2month to produce monthly tas and precipitation gridded data for all land cells."""

    # Defining the R script and loading the instance in Python
    r = robjects.r
    r['source'](r_script)

    # Loading the function we have defined in R
    r_generate_monthly_climate_fields = robjects.globalenv['generate_monthly_climate_fields']

    # generate monthly climate fields
    monthly_climate_data = r_generate_monthly_climate_fields(emulator_file,
                                                             tgav_file,
                                                             seed_value,
                                                             ngrids,
                                                             n_years,
                                                             start_yr,
                                                             scenario,
                                                             alpha_fractions_file)

    return monthly_climate_data.rx2('tas'), monthly_climate_data.rx2('pr')


def field_to_xanthos(field,
                     xth_coords_file,
                     start_yr,
                     end_yr,
                     n_years,
                     model,
                     scenario,
                     realization,
                     setting,
                     climate_data_output_dir):
    """Convert an2month downscaled climate data to the format and units required by xanthos.

    NOTE:  A single grid cell off the coast of New Zealand was off by a single grid position
    from what is represented in the climate data. This is corrected in the code.

    """

    # convert units and shape
    if setting == 'pr':
        field_arr = convert_pr(field.rx2('data'), start_yr, n_years)
        units = 'mmpermth'

    elif setting == 'tas':
        field_arr = field.rx2('data').T - 273.15
        units = 'celsius'

        # get coordinates from field object
    with localconverter(robjects.default_converter + pandas2ri.converter):
        coords_df = robjects.conversion.rpy2py(field.rx2('coordinates'))

    # get xanthos coordinates
    xth_coords_df = pd.read_csv(xth_coords_file)

    # correct new zealand grid cell that is off by a grid
    xth_coords_df['latitude'] = np.where((xth_coords_df['latitude'] == -49.75) & (xth_coords_df['longitude'] == 178.75),
                                         -49.25,
                                         xth_coords_df['latitude'])

    # merge coordinate data via left join to xanthos data
    xdf = pd.merge(xth_coords_df[['grid_id', 'latitude', 'longitude']],
                   coords_df,
                   left_on=['latitude', 'longitude'],
                   right_on=['lat', 'lon'],
                   how='left')

    assert xdf.loc[xdf.column_index.isnull()].shape[0] == 0, "NaN values found in coordinate join."

    # sort df by the xanthos grid id
    xdf = xdf.sort_values(by='grid_id')

    # reorder field data based on Xanthos grid order
    new_field = field_arr[xdf.column_index.values - 1, :]

    if climate_data_output_dir is not None:
        file_name = f'monthly_{setting}_{units}_{model}_{scenario}_{start_yr}_{end_yr}_{realization}.npy'
        output_path = os.path.join(climate_data_output_dir, file_name)
        print(f"Writing monthly climate file:  {output_path}")
        np.save(output_path, new_field)

    return new_field


def generate_drought_statistics(xanthos_config_file,
                                xanthos_output_dir,
                                xanthos_drought_thresholds_file,
                                xanthos_pet,
                                xanthos_runoff,
                                model,
                                scenario,
                                realization,
                                pr_field,
                                tas_field,
                                start_yr,
                                end_yr,
                                output_gridded_runoff=True,
                                output_runoff_by_basin=True):
    """Run xanthos to generate drought statistics for severity, intensity, and frequency."""

    # initialize config file
    xth = Xanthos(xanthos_config_file)

    # additional arguments for xanthos that override config file settings
    args = {}
    args['PrecipitationFile'] = pr_field
    args['trn_tas'] = tas_field
    args['ProjectName'] = f"{xanthos_pet}_{xanthos_runoff}_{model}_{scenario}_{realization}"
    args['OutputNameStr'] = f"{xanthos_pet}_{xanthos_runoff}_{model}_{scenario}_{realization}"
    args['OutputFolder'] = os.path.join(xanthos_output_dir, f'{xanthos_pet}_abcd_{model}_{scenario}')
    args['drought_thresholds'] = xanthos_drought_thresholds_file
    args['StartYear'] = start_yr
    args['EndYear'] = end_yr

    # create a file for gridded runoff
    if output_gridded_runoff:
        args['output_vars'] = 'q'

    # create a file that aggregates runoff by basin
    if output_runoff_by_basin:
        args['AggregateRunoffBasin'] = 1
    else:
        args['AggregateRunoffBasin'] = 0

    # run xanthos
    xth_results = xth.execute(args)


def run_drought_workflow(realization,
                         model,
                         scenario,
                         n_years,
                         start_yr,
                         end_yr,
                         emulator_file,
                         tgav_file,
                         seed_value,
                         n_fields,
                         alpha_fractions,
                         r_script,
                         xth_coords_file,
                         xanthos_config_file,
                         xanthos_output_dir,
                         xanthos_drought_thresholds_file,
                         xanthos_pet,
                         xanthos_runoff,
                         output_gridded_runoff=True,
                         output_runoff_by_basin=True,
                         climate_data_output_dir=None):
    """Run the full workflow to calculate drought statistics for severity, intensity, and frequency
    from xanthos using climate fields generated by fldgen and downscaled by an2month.

    """

    # run fldgen to generate yearly climate fields and then downscale with an2month
    tas_field, pr_field = generate_monthly_climate_fields(emulator_file=emulator_file,
                                                          tgav_file=tgav_file,
                                                          seed_value=seed_value,
                                                          ngrids=n_fields,
                                                          n_years=n_years,
                                                          start_yr=start_yr,
                                                          scenario=scenario,
                                                          alpha_fractions_file=alpha_fractions,
                                                          r_script=r_script)

    # convert to the format expected by xanthos
    pr_field = field_to_xanthos(field=pr_field,
                                xth_coords_file=xth_coords_file,
                                start_yr=start_yr,
                                end_yr=end_yr,
                                n_years=n_years,
                                model=model,
                                scenario=scenario,
                                realization=realization,
                                setting='pr',
                                climate_data_output_dir=climate_data_output_dir)

    tas_field = field_to_xanthos(field=tas_field,
                                 xth_coords_file=xth_coords_file,
                                 start_yr=start_yr,
                                 end_yr=end_yr,
                                 n_years=n_years,
                                 model=model,
                                 scenario=scenario,
                                 realization=realization,
                                 setting='tas',
                                 climate_data_output_dir=climate_data_output_dir)

    # run xanthos to generate drought statistics
    generate_drought_statistics(xanthos_config_file=xanthos_config_file,
                                xanthos_output_dir=xanthos_output_dir,
                                xanthos_drought_thresholds_file=xanthos_drought_thresholds_file,
                                xanthos_pet=xanthos_pet,
                                xanthos_runoff=xanthos_runoff,
                                model=model,
                                scenario=scenario,
                                realization=realization,
                                pr_field=pr_field,
                                tas_field=tas_field,
                                start_yr=start_yr,
                                end_yr=end_yr,
                                output_gridded_runoff=output_gridded_runoff,
                                output_runoff_by_basin=output_runoff_by_basin)


if __name__ == "__main__":

    args = sys.argv

    # SLURM_ARRAY_TASK_ID
    slurm_idx = int(args[1])

    # project level settings
    model = args[2]  # 'GFDL-ESM2M'
    scenario = args[3]  # 'rcp26'

    # if you want to run 10,000 realizations of GCM-RCP, and want to submit
    #  a job running 100 realizations per node using 100 nodes via a SLURM_ARRAY_TASK_ID
    #  you would set 'n_realizations' to 10000 and 'batch_size' to 100 and use
    #  '--array=0-99' to spread over 100 nodes.
    n_realizations = int(args[4])
    batch_size = int(args[5])

    # create a list of lists containing the number of iterations per node
    subset_realizations = get_realizations(n_realizations, batch_size, slurm_idx)

    start_yr = 1861
    end_yr = 2099
    n_years = end_yr - start_yr + 1

    # fldgen inputs
    emulator_file = f'/rcfs/scratch/d3y010/fldgen-{model}.rds'
    tgav_file = f'/rcfs/scratch/d3y010/fldgen-{model}_{scenario}.csv'

    # an2month alpha fractions
    alpha_fractions = f'/rcfs/scratch/d3y010/alpha_{model}_{scenario}.rds'

    # output directory for monthly downscaled climate data if desired
    climate_data_output_dir = None

    # R code to run fldgen and an2month downscaling
    r_script = '/people/d3y010/projects/lemens/code/fldgen_an2month.R'

    # template configuration file for xanthos
    xanthos_config_file = '/people/d3y010/projects/lemens/xanthos_lemens_template.ini'

    # coordinate reference file for xanthos
    xth_coords_file = '/rcfs/scratch/d3y010/xanthos_0p5deg_landcell_reference.csv'

    # xanthos output directory
    xanthos_output_dir = '/people/d3y010/projects/lemens/outputs/xanthos'

    # directory containing xanthos drought thresholds
    xanthos_drought_thresholds_file = f'/rcfs/scratch/d3y010/drought_thresholds_{model}_16610101-20991231.npy'

    # xanthos PET module
    xanthos_pet = 'trn'
    xanthos_runoff = 'abcd'

    print("Using the following setup:")
    print(f"model:  {model}")
    print(f"scenario:  {scenario}")
    print(f"start_yr:  {start_yr}")
    print(f"end_yr:  {end_yr}")
    print(f"n_years:  {n_years}")
    print(f"emulator_file:  {emulator_file}")
    print(f"tgav_file:  {tgav_file}")
    print(f"alpha_fractions:  {alpha_fractions}")
    print(f"climate_data_output_dir:  {climate_data_output_dir}")
    print(f"r_script:  {r_script}")
    print(f"xanthos_config_file:  {xanthos_config_file}")
    print(f"xth_coords_file:  {xth_coords_file}")
    print(f"xanthos_output_dir:  {xanthos_output_dir}")
    print(f"xanthos_drought_thresholds_file:  {xanthos_drought_thresholds_file}")
    print(f"xanthos_pet:  {xanthos_pet}")
    print(f"xanthos_runoff:  {xanthos_runoff}")

    # run the workflow
    for realization in subset_realizations:

        # generate random seed for the climate fields
        seed_value = generate_random_integer_seed()

        print(f"Processing realization:  {realization} with seed value {seed_value}.")

        run_drought_workflow(realization=realization,
                             model=model,
                             scenario=scenario,
                             n_years=n_years,
                             start_yr=start_yr,
                             end_yr=end_yr,
                             emulator_file=emulator_file,
                             tgav_file=tgav_file,
                             seed_value=seed_value,
                             n_fields=1,
                             alpha_fractions=alpha_fractions,
                             r_script=r_script,
                             xth_coords_file=xth_coords_file,
                             xanthos_config_file=xanthos_config_file,
                             xanthos_output_dir=xanthos_output_dir,
                             xanthos_drought_thresholds_file=xanthos_drought_thresholds_file,
                             xanthos_pet=xanthos_pet,
                             xanthos_runoff=xanthos_runoff,
                             output_gridded_runoff=True,
                             output_runoff_by_basin=True,
                             climate_data_output_dir=climate_data_output_dir)
        collect_garbage()
