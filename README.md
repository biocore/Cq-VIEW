# Cq-VIEW
A qPCR-based companion to the [C-VIEW](https://github.com/ucsd-ccbb/C-VIEW) (COVID-19 VIral Epidemiology Workflow) 
pipeline that calculates samples' viral loads. 

# Table of Contents
1. [Overview](#overview)
2. [Installing the Software](#installing-the-software)
3. [Configuring Cq-VIEW](#configuring-cq-view)
4. [Running Cq-VIEW](#running-cq-view)

## Overview

Cq-VIEW, a qPCR-based companion to the [C-VIEW](https://github.com/ucsd-ccbb/C-VIEW) (COVID-19 VIral Epidemiology Workflow) 
pipeline, calculates samples' viral loads.  The software includes specialized handling to support the
[SEARCH dashboard](https://searchcovid.info/dashboards/wastewater-surveillance/) as well as the UCSD campus dashboard and several others.  Cq-VIEW also outputs a standardized data file format for all target quantitation 
methods (such as qPCR and ddPCR) so that their data can be used interchangeably
in downstream analyses.

Full intermediate and results files are stored locally with the 
following folder-/file-naming structure.  Throughout, "sample" refers
to a particular qPCR measurement of a particular extraction of a particular physical specimen 
(e.g., `3.24.23.PLMAR21.R1.A15`) while "sample_base" refers to a particular site on
a particular day (e.g., `PLMAR21`)

* <local_reports> folder (e.g., `/Users/me/covid_wastewater_qpcr/reports`)
    * <report_name> folder (e.g., `2023-04-20_12-13-10_230327_RTL`), made up of analysis timestamp followed by report_base_name (usually the qPCR run name)
      * <report_run_name> folder: either <report_name>_only (e.g., `2023-04-20_12-13-10_230327_RTL_only`) or <report_name>\_cumulative\_since\_\<date> (e.g., `2023-04-20_12-13-10_230327_RTL_cumulative_since_2023-03-16`)
        * `generate_reports_error.log`: overall error log file; this should always be empty, and report results should not be used if it is not
        * `input` folder:
          * One or more ThermoFisher QuantStudio 7 Pro qPCR output files in xlsx or csv format (e.g., `230327_RTL_Results_20230328 125315.csv`)
          * [Optional] One or more `exclude` files (e.g., `230327_RTL_exclude.yml`) specifying data from a particular run to be excluded from analysis
        * `output` folder
          * One folder per `site group` (e.g., `ucsd`, `usd`, `sandiego`, `ucsf`, etc)
            * verbose logs for each step of the analysis
            * <report_run_name>\_<site_group>\_ladder\_\<target>.csv (e.g., `2023-04-20_12-13-10_230327_RTL_only_ladder_Promega_N1.csv`): file of ladder values for each target
            * <report_run_name>\_<site_group>\_ladder\_curves\_\<assay>.pdf (e.g., `2023-04-20_12-13-10_230327_RTL_only_ladder_curves_Promega_SARS-CoV-2_N1`): plots of the curves derived from each target ladder, along with the standard curve for that target (as found in the `config.yml`)
            * [Only for `sandiego`] `input` folder:
              * `*_sewage_seqs.csv` files: current files downloaded from the [SEARCH dashboard repo](http://www.github.com/andersen-lab/SARS-CoV-2_WasteWater_San-Diego)
            * `output` folder:
              * <report_run_name>_<site_group>.xlsx (e.g., `2023-04-20_12-13-10_230327_RTL_only_ucsd.xlsx`): human-readable spreadsheet containing viral load calculations for each sample and each sample base
              * <report_run_name>\_<site_group>\_sample_data.csv: calculations of viral load for each sample
              * <report_run_name>\_<site_group>\_sample_base_data[_\<site>].csv (e.g., `2023-04-20_12-13-10_230327_RTL_only_sandiego_sample_base_data_PL.csv`): calculations of viral load for each sample base, optionally limited to a particular site within the site group (e.g., the `PL` site within the `sandiego` site group)
              * <report_run_name>\_<site_group>\_dashboard_data.csv: per-sample viral load information formatted for delivery to the CMI dashboard
              * <report_run_name>\_<site_group>\_standard_data.csv: per-sample viral load information formatted for interoperability with data from other quantitation approaches.
              * [Only for `sandiego`]:
                * `*_concs.pdf` and `*_cqs.pdf` files (e.g., `2023-04-20_12-13-10_230327_RTL_only_sandiego_concs.pdf`):  plots visualizing Cq values and concentrations (i.e., viral loads) of each measurement
                * `*_sewage_seqs.csv` files: Updated files containing the new viral loads from the latest run
                * <report_run_name>\_<site_group>\_repo\_upload.sh (e.g., `2023-04-20_12-13-10_230327_RTL_only_sandiego_repo_upload.sh`): generated script that, when run, creates a pull request to submit the updated `*_sewage_seqs.csv` files to the [SEARCH dashboard repo](http://www.github.com/andersen-lab/SARS-CoV-2_WasteWater_San-Diego)

## Installing the Software

### SEARCH Alliance prerequisites (optional)

Using Cq-VIEW to process SEARCH Alliance data requires the following information:
* A Github account with known user name and email,
  * Having contributor access to https://github.com/AmandaBirmingham/SARS-CoV-2_WasteWater_San-Diego
  * Attached to an ssh key, with known public and private keys

Install SEARCH-associated prerequisites:
   1. Configure `git` check-in information
      1. Run `git config --global --edit` to edit the config doc; then
      2. Set the user name and email for the chosen github account, uncomment those lines, and resave
   2. Configure Github ssh access
      1. Get the RSA key-pair (both public and private) for github user's ssh key
      2. Copy the key-pair into the `.ssh` folder in the user's home directory
      3. Ensure the private key has limited permissions by running `chmod 600 <keyname>`
      4. Validate ssh access to github is working
         1. Run `ssh -T git@github.com`
         2. Confirm the addition of github.com to the list of known hosts
         3. Access is validated if you receive a message stating `You've successfully authenticated, but GitHub does not provide shell access`
   3. Clone the SEARCH repo (note that this must come AFTER configuring Github ssh access)
   
      1. Run `git clone git@github.com:AmandaBirmingham/SARS-CoV-2_WasteWater_San-Diego.git` to clone the (public) fork of the Andersen lab repo for dashboard inputs
      3. Run `cd SARS-CoV-2_WasteWater_San-Diego`
      4. Run `git remote add andersen_wastewater https://github.com/andersen-lab/SARS-CoV-2_WasteWater_San-Diego.git` to add the Andersen lab repo as an upstream to the fork repo

### Cq-VIEW

Install Cq-VIEW and its dependencies:
   4. Install the `wastewater` conda environment for Cq-VIEW
   
      1. Move out of the `SARS-CoV-2_WasteWater_San-Diego` directory (if necessary)
      2. Download the `wastewater_env.yml` file from this repo
      3. Run `conda env create -f wastewater_env.yml`
      4. Run `conda activate wastewater`
   5. Install the `xlsxwriter` package into the `wastewater` conda environment using `pip install xlsxwriter`
   6. Install Cq-VIEW from github
      
      1. Run `git clone https://github.com/biocore/Cq-VIEW.git`
      2. Run `pip install -e .`

If processing SEARCH Alliance data, also add `gh` to new conda enviroment:
   7. Configure `cli` access for `gh`
   
      1. Run `gh auth login`
      3. Choose `GitHub.com` as the account to log into
      4. Choose `SSH` as the protocol
      5. Choose the ssh key installed above
      6. Title the key `<keyname>`
      7. Open a browser and *ensure that you are logged in with your chosen github account!*
      8. As directed by the `gh` auth process, open the url `https://github.com/login/device`
      9. Paste in the one-time code provided by the process and confirm access
         1. If you receive an "[HTTP 422 error key is already in use](https://github.com/cli/cli/discussions/6271)" error message, this is ok, as described in the link

## Configuring Cq-VIEW

Much of Cq-VIEW's behavior is specified by the configuration yaml file (`config.yml` for all processing except that for UCSF, which uses `config_ucsf.yml`).  

### Selecting the Curve Model

Cq-VIEW calculates viral loads from Cq values using an equation derived from one or more standard curves.  The software supports three different kinds of curve models: single, pooled, and mixed.  These model types are defined in ["Improved strategies and optimization of calibration models for real-time PCR absolute quantification"](https://pubmed.ncbi.nlm.nih.gov/20701947/) (Sivaganesan et al, *Water Res.* 2010, PMID: 20701947).  In short:

* The "single" model derives a completely new curve for every run, using ladders from that run.  
* The pooled model derives one curve from a pre-existing set of ladders and uses that curve for all future runs.  
* The mixed model derives the slope of the  curve from a pre-existing set of ladders but derives the intercept for each run from a single dilution point in the ladder included on that run.  

To select the curve model used by Cq-VIEW, find the `analysis_settings` section of the config file and set the value of its `curve_model` key to be `single`, `pooled`, or `mixed`:

```
"analysis_settings":
  "curve_model": pooled  # options: single, pooled, mixed
```

If the value is `pooled` or `mixed`, it is also necessary to specify the pre-derived standard curve for each target in the `standard_curves` section of the config, as in the following example:

```
"standard_curves":
  "PMMoV_N1":
    "curve_name": "2023-02-15_15-25-02_1_curves"
    "slope": -3.500070
    "intercept": 36.417650
    "use_ln": False
    "reaction_vol_per_l": 1000000
```

In addition, the `ladders` part of the `input_files` section contains further settings defining the lower bound of the ladder included in each run file, as well as the number of dilution steps done, as shown here:

```
"input_files":
  "ladders":
    "ladder_bottom": 0.004
    "ladder_steps": 8
    "intercept_ladder_step": 3
    "regression_bottom_step": 3
```

Only steps at or above the `regression_bottom_step` (when counting from the bottom) will be included when calculating the standard curve from a run (since very low dilutions tend to fall off the linear part of the curve). Additionally, the `intercept_ladder_step` setting is used by the `mixed` model to specify the single dilution point from run's ladder to use when deriving the intercept for that run.   
      
## Running Cq-VIEW

*Note: uploading to any campus dashboard requires an AWS user with S3 
permissions for that dashboard, with a known access key and secret access key.*

1. Activate the `wastewater` conda environment
2. For qPCR data:
   1. Load the qPCR result file(s)
      1. Place each qPCR results file provided by the lab in a new, clearly named run directory
         1. For example, put `230327_RTL_Results_20230328 125315.csv` into `/Users/me/covid_wastewater_qpcr/data/230327_RTL`
      2. Capture the name(s) of the run director(ies) (e.g., `230327_RTL`) for use below
   2. For UCSD or San Diego data (but **not** UCSF data), run `generate_report` 
      1. Input the following arguments:
         1. <conda_bin_dir> (e.g., `"/Applications/miniconda3/bin"`): the `bin` directory for the local `conda`/`Anaconda` installation
         2. <search_repo_parent_dir> (e.g., `"/Users/me/Repositories"`): the local folder containing the SEARCH Alliance `SARS-CoV-2_WasteWater_San-Diego` repository clone
         3. <qpcr_data_parent_dir> (e.g., `"/Users/me/covid_wastewater_qpcr/data"`): the local folder containing the new qPCR results directory created above
         4. <report_parent_dir> (e.g.,`"/Users/me/covid_wastewater_qpcr/reports"` ): the local folder into which the new report output folder will be created
         5. <report_base_name> (e.g., "`test_report`; may be `""`): a custom name for the report; if an empty string is given, then a concatenation of the source directory names is used
         6. <earliest_date> (e.g., `"05/01/2023"`; may be `""`): the earliest date for which to include existing qPCR data in the analysis; ignored for non-`sandiego` site groups
            1. Pick a date that includes the dates sampled in the next-most recent batches of `sandiego` samples (e.g. usually a week ago or so) unless you specifically want to see a comparison of the results over a longer time-period
         7. <run_directories_names> (e.g., `"230512_RTL_USD_PL_ENC"`; may include multiples separated by spaces, like `"230327_RTL" "230512_RTL_USD_PL_ENC"`): the name(s) of the directories within <report_parent_dir> whose results should be included

         ```
         # Example:
      
         conda activate wastewater
         generate_report "/Applications/miniconda3/bin" \
           "/Users/me/Repositories" \
           "/Users/me/covid_wastewater_qpcr/data" \
           "/Users/me/covid_wastewater_qpcr/reports" \
           "" \
           "05/01/2023" \
           "230512_RTL_USD_PL_ENC" 
         ```
      2. Note that if UCSF data is run through the `generate_report` function, it will fail with the error `ValueError: Should not process ucsf data with default config!`
         1. This is a safeguard to ensure that no UCSF data is processed with the wrong settings, because UCSF data must be processed with a different config file than non-UCSF data
         2. If you receive this error, you are trying to do something illegal; instead, follow the instructions for processing UCSF data
   3. For UCSF data (but **not** UCSD or San Diego data), run `generate_ucsf_report`
      1. Input the following arguments:
         1. <conda_bin_dir> (e.g., `"/Applications/miniconda3/bin"`): the `bin` directory for the local `conda`/`Anaconda` installation
         2. <search_repo_parent_dir> (e.g., `"/Users/me/Repositories"`): the local folder containing the SEARCH Alliance `SARS-CoV-2_WasteWater_San-Diego` repository clone
         3. <qpcr_data_parent_dir> (e.g., `"/Users/me/covid_wastewater_qpcr/data"`): the local folder containing the new qPCR results directory created above
         4. <report_parent_dir> (e.g.,`"/Users/me/covid_wastewater_qpcr/reports"` ): the local folder into which the new report output folder will be created
         5. <report_base_name> (e.g., "`test_report`; may be `""`): a custom name for the report; if an empty string is given, then a concatenation of the source directory names is used
         7. <run_directories_names> (e.g., `"230607_UCSF"`; may include multiples separated by spaces, like `"230607_UCSF" "230609_UCSF"`): the name(s) of the directories within <report_parent_dir> whose results should be included

         ```
         # Example:
      
         conda activate wastewater
         generate_ucsf_report "/Applications/miniconda3/bin" \
           "/Users/me/Repositories" \
           "/Users/me/covid_wastewater_qpcr/data" \
           "/Users/me/covid_wastewater_qpcr/reports" \
           "" \
           "230607_UCSF" 
         ```
      2. Note that if non-UCSF data is run through the `generate_ucsf_report` function, it will fail with the error `ValueError("Should not process non-ucsf data with non-default config!")`
         1. This is a safeguard to ensure that no non-UCSF data is processed with the wrong settings, because UCSF data is processed with a different config file than non-UCSF data
         2. If you receive this error, you are trying to do something illegal; instead, follow the instructions for processing non-UCSF data
   4. Evaluate and manually deliver results
      1. Open the <report_name>_only folder 
      2. Ensure that the `generate_reports_error.log` is 0 B
      3. Open the `output` folder
      4. Open the <report_run_name>\_<site_group>\_ladder\_curves\_\<assay>.pdf file and examine the ladder plots to ensure that neither of them has a radically different slope than the standard curve associated with that target
      5. Handle each sub-folder in `output` based on its site-group:
         1. `usd`:
            1. Open the usd/output/<report_run_name>\_\<site_group>.xlsx file and skim the contents of the first sheet (`sb_<etc>`), which contains sample base info.  Lots of -1s and a few large positive numbers for `normalized_concentration` are usual.  If you see unusually small numbers of sample bases (<10) or more than a handful of `TRUE` values in the `ignore` column, something may be wrong and further investigation is warranted.
            2. Assuming no issues have been identified, transfer the <report_run_name>\_<site_group>\_sample_base_data.csv file into the usd AWS S3 bucket
         2. `ucsd`
            1. Open the ucsd/output/<report_run_name>\_\<site_group>.xlsx file and skim the contents of the first sheet (`sb_<etc>`), which contains sample base info.  Lots of -1s and a few large positive numbers for `normalized_concentration` are usual.  If you see unusually small numbers of sample bases (<30) or more than a handful of `TRUE` values in the `ignore` column, something may be wrong and further investigation is warranted.
            2. Assuming no issues have been identified, transfer the <report_run_name>\_<site_group>\_**dashboard**_data.csv file into the ucsd AWS S3 bucket
         3. `ucsf`:
            1. Open the ucsf/output/<report_run_name>_ucsf_delivery.csv file and skim the contents and skim the `target_to_control_ratio` values.  Numbers in the 1E-3 to 1E-4 range are typical.  If you see larger or smaller values than this, or more than a couple of NaNs, further investigation may be warranted.
            3. Email the *_ucsf_delivery.csv to the UCSF contacts.
         4. `sandiego`:
            1. Open the usd/output/<report_run_name>\_\<site_group>.xlsx file and skim the contents of the first sheet (`sb_<etc>`), which contains sample base info. There are usually only a few sample bases (2-8), and none are expected to have -1 for `normalized_concentration`.  If they do, something may be wrong and further investigation is warranted.
         
               1. For any sample base with `ignore = TRUE` that does NOT have the QC message `Is not within 2-fold of a context value for normalized_concentration`, check on the sample tab (`s_<etc>`) to see what issues affected that sample base's constituent samples.  In most such cases, the lab is asked to rerun the relevant sample(s) (one time only) to see if the sample base can be recovered.
               2. Any sample base with `ignore = TRUE` that DOES have the QC message `Is not within 2-fold of a context value for normalized_concentration` does not receive further examination; it may be "reclaimed" in the cumulative report or in a subsequent run if future sample bases from the same site have `normalized_concentration` values within range of it
            3. Note which sample bases are included in this qPCR run, for use later
      6. If SEARCH Alliance data is present (i.e. if and only if there is a `sandiego` folder), go up to the top level of the report folders and open the <report_name>\_cumulative\_since\_\<date> folder (this is only produced if `sandiego` samples are in the qPCR run)
         1. Ensure that the `generate_reports_error.log` is 0 B
         2. Open the `output/sandiego/output` folder
         3. Open the `*_concs.pdf` and `*_cqs.pdf` plot files and examine the Cqs and concentrations of the latest `sandiego` samples in context of other recent samples from the same sites within the site group.  If they look notably different, stop and discuss with the lab.
         4. Run the <report_run_name>\_<site_group>\_repo\_upload.sh bash script to commit the updated `*_sewage_seqs.csv` files to the [https://github.com/AmandaBirmingham/SARS-CoV-2_WasteWater_San-Diego](https://github.com/AmandaBirmingham/SARS-CoV-2_WasteWater_San-Diego) public fork of the SEARCH dashboard repo and then make a pull request for them in the actual [SEARCH dashboard repo](http://www.github.com/andersen-lab/SARS-CoV-2_WasteWater_San-Diego)
         5. Examine the new pull request in [SEARCH dashboard repo](http://www.github.com/andersen-lab/SARS-CoV-2_WasteWater_San-Diego) and check:
            1. That no previously reported values have been changed (this is a rare unintended side effect of including past data using the `earliest_date` approach)
            2. Whether any previously excluded sample bases from previous runs have been added due to now passing context QC (this is intended behavior, but good to notice)
            3. That the added values look sensible
         6. Contact Josh Levy in the Andersen lab and inform him a PR is ready for review for the SEARCH dashboard repo
3.  For ddPCR data:
    1. Ensure the ddPCR result file(s) contain the required per-result metadata in the built-in columns
       1. `Sample description 1`: unique sample name (e.g., `PL_8_17_22_3`)
       2. `Sample description 2`: name of sampler used to collect sample (e.g., `Point Loma`)
       3. `Sample description 3`: sample collection date (e.g. `8/17/22`)
       4. `Sample description 4`: unique name of this quantitation result on this sample (e.g., `PL_8_17_22_2.F08.1021083101G21000050`)
    2. Run `generate_ddpcr_report`
       1. Input the following arguments:
          1. <report_parent_dir> (e.g.,`"/Users/me/covid_wastewater_qpcr/reports"` ): the local folder into which the new report output will be created
          2. <protocol_name> (e.g., `"ddPCR with ddPCR Multiplex Supermix"`): the name of the quantification protocol used to generate the data
          3. <protocol_date> (e.g., `"2022-09-13"`): the date on which the quantification protocol was performed
          4. <control_target_name> (e.g.,`"PMMoV"` ): the name of the target to be used as a fecal control
          5. <report_base_name> (e.g., "`test_report`; may be `""`): a custom name for the report; if an empty string is given, then a concatenation of the source file names is used
          6. <quant_file_paths> (e.g., `"/Users/me/Desktop/220913_samples_data.csv"`; may include multiples separated by spaces, like `"/Users/me/Desktop/220913_samples_data.csv" "/Users/me/Desktop/221005_samples_data.csv"`): the path(s) to the files whose results should be included

      ```
      # Example:
      
      conda activate wastewater
      generate_ddpcr_report \
        "/Users/me/covid_wastewater_ddpcr/reports" \
        "ddPCR with ddPCR Multiplex Supermix" \
        "2023-09-13" \
        "PMMoV" \
        "" \
        "/Users/me/Desktop/220913_samples_data.csv" 
    ```
    3. Evaluate results:
       1. Look for the `*_standard_data.csv` in the specified report parent directory
       2. Note the `external_target_name` column, since multiple targets can be evaluated for the same sample
       2. Examine `target_copies_per_liter`, which is the primary ddPCR output
          1. This will be a number greater than or equal to zero if the target was successfully quantified
          2. If the target was NOT successfully evaluated, this column will be empty and the `accept` column will be `False`
       3. Examine `target_to_control_conc_ratio`
          1. This will be a number if both the target in question and the fecal control were successfully quantified for this sample
          2. If either was NOT successfully evaluated, this column will be empty
       2. Note that the standard `target_cq` column will always be empty for ddPCR

