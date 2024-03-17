# AutomatedTargetScreening
Source code, file descriptions, and installation instructions for the AutomatedTargetScreening (ATS) workflow. 

---

# Installation instruction
The package is based on R version 4.1.3. The package is installed using the `devtools` R package:

```R
devtools::install_github("braungeorg/AutomatedTargetScreening")
```

This will install ATS and its function documentation as well as CRAN-accessible dependencies to your local R and R Studio. However, the R package `rawrr` needs to be installed seperately using the R package `BiocManager`. All dependencies including `rawrr` are installed if the `Automated_Target_Screening_workflow.R` script within the example environment from Zenodo (10.5281/zenodo.10047377). However, you can also run the following code to install all packages including `rawrr`. 

```R
packages <- c("dplyr", "rawrr", "RaMS", "openxlsx","qpdf","data.table",
              "seewave","doSNOW","doParallel","lubridate","zoo","gtools","ggplot2",
              "stringr","ggtext","yaml","foreach","MESS","tcpl","jmotif","plyr","outliers")

# Function to check if a package is installed
is_installed <- function(pkg) {
  return(requireNamespace(pkg, quietly = TRUE) || pkg %in% installed.packages()[, "Package"])
}

# Check if packages are already installed
missing_packages <- packages[!sapply(packages, is_installed)]

# Install missing packages
if (length(missing_packages) > 0) {
  # Install packages from CRAN
  cran_packages <- missing_packages[!missing_packages %in% c("rawrr")]
  if (length(cran_packages) > 0) {
    install.packages(cran_packages, dependencies = TRUE)
    lapply(cran_packages, requireNamespace, quietly = TRUE)
  }
  
  # Install packages from Bioconductor
  if ("rawrr" %in% missing_packages) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("rawrr")
    library(rawrr)
    rawrr::installRawrrExe()
    rawrr::installRawFileReaderDLLs()
  }
}

#Load packages
invisible(lapply(packages, library, character.only = TRUE))
```
Please note that the functions `installRawFileReaderDLLs()` and `installRawrrExe()` of the `rawrr` package need to be executed before the ATS script can be used for `.raw` files. 

Once the package and all its dependencies are installed to your library you can use the `Environment` folder provided on Zenodo via 10.5281/zenodo.10047377. Copy this `Environment` folder to wherever you like. 

It already includes the .R scripts within the `0_Scripts` folder which can be used for analysis of sample batches included in the `2_Samples` folder. The `1_Targets` folder includes the target lists, assignments of internal standards, lists of multiple peaks, and list of analytes for qualitative analysis. There are already folders created in `2_Samples` for *LC* and *GC* and their respective .yaml files for analysis settings. These example files cover only a small set of analytes for LC-ESIpos and GC-EI. The sample files were minified .mzML files which were created using the `RaMS`R package. Further, the `3_Results` folder already includes the results after analysis with ATS. 

---

# Description of data files

## 0_Scripts

### Automated_Target_Screening_Functions

This .R file is including all functions and will be sourced by the `workflow` script. 

### Automated_Target_Screening_workflow

This .R file is the main script and can be ran individually or as batch. It is dependent on a correctly setup `environment`. 

## 1_Targets

Targets included all files describing the analytes and internal standards. 

### target_list

The target_list includes information about all quantified analytes as well as internal standards. 

- ID: Unique per compound and has different hierarchies
 - Numeric: Quantified ions, i.e. 1, 002
 - Numeric + letter: Confirming ions, i.e. 1a,002c
 - Starts with "IS": Internal standards, i.e. IS01, IS508
- m/z: The mass to charge ratio of the respective ion
- retention time: The expected retention times in min. Note that they have to be correct relatively to each other for shift calculations to be meaningful
- identitity: A combination of ID and chemical identifier (i.e. 0016_Diuron)
- Formula: Not necessary but a column which can be used for the molecular formula of the analytes
- comment: A column for comments
- ion: A way to define the ion like [M+H]+
- prim.sec: Which MS level should be used for this ion (MS1 or MS2)

### qual_list

A list of analytes which are only screened in a qualitative approach (i.e. are only included as reference standards without quantification or at a single concentration)

- ID: Unique per compound and has different hierarchies. For qualitative lists it is recommended to add another prefix like "Q" to allow separate ID pools
 - Numeric: Quantified ions, i.e. Q1, Q002
 - Numeric + letter: Confirming ions, i.e. Q1a,Q002c
 - Starts with "IS": Internal standards, i.e. IS01, IS508
- m/z: The mass to charge ratio of the respective ion
- retention time: The expected retention times in min. Note that they have to be correct relatively to each other for shift calculations to be meaningful
- identitity: A combination of ID and chemical identifier (i.e. Q016_Diuron)
- Formula: Not necessary but a column which can be used for the molecular formula of the analytes
- comment: A column for comments
- ion: A way to define the ion like [M+H]+
- prim.sec: Which MS level should be used for this ion (MS1 or MS2)

### IS_assignment

A list of the manually pre-set assignment of internal standards to analytes. 

- Compound: The identification of the compound. Needs to be equal to the `identity` of target and qualitative lists.
- ISTD: The `identity` of internal standards also used in the target lists.

### Multiple_peaks

A file defining if multiple peaks occur within a likely retention time window. Defines which peaks are finally selected for deriving a final retention time of the analyte.

- ID: Numeric ID of the analyte (i.e. 0001 -> 1).
- identitiy: Identity of the analyte from the `target_list`
- multiple: From left to right which peak is the relevant one (1 = first, 2 = 2nd, ..., 0 = all).

## 2_Samples

This folder contains all respective measurmements split by mode (i.e. ESIpos, ESIneg, GC). Each folder contains the raw data in either *.raw* or *.mzML* file format. Further, a **.yaml** file needs to be included to define settings for the analysis:

Files in `1_Targets`:
- Targets_file: Character, the name of the `target_file`
- Qual_file: Character, the name of the `qual_list`
- IS_Assignment: Character, the name of the `IS_assignment`
- Multiple_peaks_file: Character, the name of the `Multiple_peaks`

Tags/IDs of files included in `2_Samples`:
- Reference_ID: Character, identifier used for samples which contain quantified reference standards like i.e. "CAL","QA/QC".
  > *IMPORTANT: "CAL_xxxx" is the necessary syntax for quantified samples, i.e. CAL_0p01 for 0.01, CAL_00p1 for 0.1, CAL_0001 for 1 etc.*
- Qual_ID: Character, identifier which sample(s) should be used to derive peaks for qualitative-only analytes
- Solvent_Blank: Chracter, identifier for samples which were only solvent_blanks (i.e. "MeOH" for pure methanol injection)

Settings relevant for handling the mass spectrometry data:
- FullscanMS1: Character, which filter shall be used for MS1 data (relevant for *.raw* files)
- ppm_val: Numeric, what m/z error shall be allowed in ppm
- method_time: Numeric, the duration of the measurement in minutes (needed to define plausible shift ranges)
- Scan filters: Character, a list of scan filters applied for analysis (relevant for *.raw* files)

Settings relevant for the automatic data evalution:
- alphabet_size: Numeric, the range/size of the alphabet used in symbolic aggregate approximation (SAX)
- minimum_search_window: Numeric, the minimum range in minutes that is screened for peaks in both directions of the aimed retention time
- use_shifting_search_window: Boolean to allow/deny to shift/adjust search window while generating the peaklist?
- maximum_nr_of_search_window_extensions: Numeric, the maximum number of search window extensions if the shifting search window is allowed per direction (left and right)
- zigzag_trigger_threshold: Numeric, the portion of the EIC in terms of intensity values that needs to be defined as zigzag, meaning highly fluctuating intensities
- normal_background_quantile: Numeric, the quantile of the background defining intensity values used to define the baseline
- higher_background_quantile: Numeric, the quantile of the background-defining intensity values used to define the baseline if the EIC is considered problematic, i.e. exceeding the zigzag_trigger_threshold
- minimum_background_ratio: Numeric, the minimum ratio of intensity and baseline an intensity value needs to exceed to be considered a peak
- extended_baseline_factor: Numeric, this factor is applied to the baseline and all intensities below factor*baseline are considered baseline
- maximum_nr_of_stagnant_intensity_values_peaktop: Numeric, the maximum number of intensities without significant decrease after which the start and end time of the peaktop is defined
- intensity_factor_decrease_peaktop: Numeric, the factor applied to the last peak-assigned intensity which needs to be undercut within the maximum number of intensities for a decrease to be considered significant in defining the peaktop
- intensity_factor_increase_edges: Numeric, the factor applied to the last peak-assigned intensity which if exceeded triggers a density-dependent counter to assign small increases to the peak while splitting peaks at the peak edges
- intensity_factor_decrease_edges: Numeric, the factor applied to the last peak-assigned intensity which needs to be undercut within the density-dependent maximum number of datapoints for a decrease to be considered significant in defining the peak edges
- density_factor_increase_edges: Numeric, the factor applied to the intensity density (datapoints within 0.1 minute) which if exceeded by having a constant increase or lack of decrease triggers a stop of the peak detection algorithm and defines the ends of the peak edges
- density_factor_decrease_edges: Numeric, the factor applied to the intensity density (datapoints within 0.1 minute) which if exceeded by having too many stagnant intensity values triggers a stop of the peak detection algorithm and defines the ends of the peak edges
- width_smoothing: Numeric, a factor which defines how the data density is increased when smoothing. Higher values increase data density but can lead to "oversaturation"
- width_factor_background: Numeric, a factor defining how a defined search window will be extended to define the background EIC
- sample_search_window: Numeric, the time in minutes added to the start and end retention time of a defined peak in the suspect screening (to improve visualization)
- outer_search_window_multiple_peaks: Numeric, the larger time in minutes added to the start and end retention time of a defined peak to generate the peaklist if peak patterns are expected, used for different intensity ratios to filter for invalid peaks
- inner_search_window_multiple_peaks: Numeric, the smaller time in minutes added to the start and end retention time of a defined peak to generate the peaklist if peak patterns are expected, used for different intensity ratios to filter for invalid peaks
- outer_intensity_ratio_multiple_peaks: Numeric, ratio of the maximum intensity within the larger RT range expected for peak patterns
- inner_intensity_ratio_multiple_peaks: Numeric, ratio of the maximum intensity within the larger RT range expected for peak patterns
- maximum_allowed_shift: Numeric, the maximum shift allowed in minutes
- maximum_allowed_shift_ratio: Numeric, the maximum fraction of IS-dependent shift of the allowed shift, only relevant if allowed shift exceeds maximum value
- maximum_nr_of_a: Numeric, a defines the "baseline" of a peak and needs a limit so that very flat peaks or peaks with only few high intensity values are not considered valid
- minimum_nr_of_high_intensity_letters: Numeric, similar to the limitation of a you can set a minimum of high-intensity letters (dependent on your selected alhabet size) to ensure that your selected peaks are highly deviating from the baseline
- minimum_nr_of_datapoints_per_peak: Numeric, the minimum number of data points needed to define a peak
- minimum_peak_area: Numeric, the minimum area needed to define a peak
- maximum_peak_width: Numeric, the maximum peak widht in minutes
- minimum_cutoff_intensity_factor: Numeric, this factor is used to define, based on the highest intensity of a compound, at which intensity value samples will be either split (exceeding this threshold) or merged (below this threshold)
- minimum_confirming_peak_height: Numeric, the intensity value which needs to be exceeded if a peak should be considered as confirming ion
- gen.plots: Boolean, generate PDFs as results plots?
- use.MINDIST: Boolean, shall MIDNIST be applied as quality criterion for selecting valid peaks?
- max_MINDIST: Numeric, which MINDIST value should be considered as acceptable MINDIST?
- maximum_RT_tolerance_start: Numeric, the maximum retention time tolerance for the start retention times of a peak
- maximum_RT_tolerance: Numeric, the maximum retetion time tolerance for the apex retention time of a peak
- maximum_RT_tolerance_end: Numeric, the maximum retention time tolerance for the end retention times of a peak
- minimum_datapoints_per_sample_peak: Numeric, the lowest number of datapoints needed to define a peak in the suspect screening
- minimum_qualitative_threshold: Numeric, the ratio of the total number of data points which needs to exceed the in-sample LOD to trigger at least qualitative detection
- use.area: Boolean, shall area (T) or height (F) be used for analysis?
- max_calibration: Numeric, ratio of highest accepted calculated value from nominal concentration. I.e. 1.2 = 1200 ng/mL are still valid if 1000 ng/mL was used as highest concentration.
- minimum_r2: Numeric, the lowest accepted coefficient of determination of calibration curves
- minimum_r2_intensity_dependent_shift: Numeric, the lowest accepted coefficient of determination for linear regression used to identify valid intensity-dependent shifts within the calibration
- maximum_error_of_quantitation: Numeric, the maximum accepted error of quantitation in percent
- minimum_nr_of_datapoints_calibration: Numeric, the minimum number of datapoints needed to define a valid calibration
- quantified_unit: Character, the concentration unit
- IS_deviation: Numeric, what deviation is acceptable for internal standards? I.e. 0.2 = 20% deviation

## 3_Results

Folder which gathers all results per samples included in `2_Samples`. Each results folder consists of: 

- plots: A folder which contains PDFs of analyte-sample pairs. Only relevant if a manual re-evaluation is considered.
- results_plots: PDFs per analyte containing all chromatograms for the whole sample batch. Only for the main/quantified ion.
- "for_manual_evaluation": A data file which shall be used if manual re-evaluation is considered
  - Sheet "What to do" give a brief introduction to valid actions and their functionality
  - Sheet "worksheet" is to be used for evaluation.
    - Compound: The identity of the analyte as used in the `target_list`
    - File: The name of the data file
    - Potential_value: Values which were calculated for this analyte-sample pair
    - Current_value: Either numeric if no problem occurred or character if i.e. "Masked_by_Background" or "CHECK"
    - Action: As defined in the sheet "What to do"
    - Start_RT: If manually set to integration, the start of the peak in min
    - End_RT: If manually set to integration, the end of the peak in min
- "Intensity_dependent_shifts": A file listing all compounds for which an intensity-dependent shift was identified. Each analyte is listed in a separate sheet.
   - int: The intensity values
   - shift: The shift in minutes from the expected retention time
- "Results_Algorithm": The main results file containing several sheets of information.
   - A data sheet for Areas/Heights (Areas/Heigths), Response Ratios (Ratios), quantified values (Final), filtered quantified values (Final_excluded), retention times
      - Compound: Character, the identity of the analyte
      - Comment: Character, comments created by ATS, mainly if ions were not found or why samples where excluded
      - mz: Numeric, the main m/z of the analyte
      - RT_start: Numeric, the peak start in minutes
      - RT: Numeric, the apex retention time in minutes
      - RT_end: Numeric, the peak end in minutes
      - ISTD: Character, the assigned internal standard
      - CV IS: Numeric, the coefficient of variance of the internal standard over the sample batch
      - Valid references: Numeric, only relevant if one-level references (i.e. for recovery) are used. Defined the number of valid references used in the analysis.
      - Valid calibration: Boolean, used to filter valid/invalid calibrations.
      - R2_calbiration: Numeric, the coefficient of determination of the respective calibration
      - Min_calibration: Numeric, the lowest nominal concentration value included in the respective calibration
      - Max_calibration: Numeric, the highest nominal concentration value included in the respective calibration
      - Matrix correction: Numeric, the ratio of matrix sample to reference sample and effect correction applied to samples
      - Bank values: Numeric, the respective mean value (area/height, response ratio, quantified) found in solvent blanks
      - Names of the samples: Numeric, the respective value (area/heigth, response ratio, quantified)
  - "Qualitative", a sheet used to summarize the occurrance of analytes in TRUE/FALSE
      - Compound: Character, the identity of the analyte
      - Names of the samples: Boolean, TRUE/FALSE if analyte is present or not (Background if masked by background and hence inconclusive)
  - "Peaklist", a summary and details of identified peaks
      - Compound: Character, the identity of the analyte
      - mz: Numeric, the main m/z of the analyte
      - Comment: Character, comments created by ATS, manly if ions were not found or why samples where excluded
      - RT_start: Numeric, the peak start in minutes
      - RT: Numeric, the apex retention time in minutes
      - RT_end: Numeric, the peak end in minutes
      - Start_RT_level: Numeric, the retention time in minutes after peak is exceeding the relative intensity level (i.e. 50% maximum intensity)
      - End_RT_level: Numeric, the retention time in minutes after peak is below the relative intensity level (i.e. 50% maximum intensity)
      - Sequence: Character, the SAX sequence of the peak
      - Nr_of_Points: Numeric, the number of data points defining the peak without smoothing
      - Width: Numeric, the peak width in minutes at a respective intensity level (i.e. 50% of maximum intensity)
      - Height: Numeric, the heigth/intensity of the peak
      - Area: Numeric, the integrated area of the peak
      - level: Numeric, the relative intensity value used to define the core peak (i.e. 0.5 = 50% maximum intensity)
      - Concentration: Numeric, the concentration level of the highest reference compound
      - Bg_Start: Numeric, the start in minutes of the the extracted ion chromatogram used for background and baseline calculations
      - Bg_End: Numeric, the end in minutes of the extracted ion chromatogram used for background and baseline calculations
      - density: Numeric, the number of smoothed data points defining 0.1 minute of the highest reference standard
      - Start_var: Numeric, accepted variability of the starting retention times (i.e. 0.3 = 30%)
      - End_var: Numeric, accepted variability of the ending retention time (i.e. 0.3 = 30%)
      - max_MINDIST: Numeric, a threshold calculated as maximum accepted MINDIST
      - RT_tol: Numeric, accepted variability of the apex retention time (i.e. 0.2 = 20%)
- "Calibration_Curves": Only for calibrations. A PDF visualizing the final calibration curves per analyte. Excluded values are red and the coefficient of determination is also displayed.
- "Reference_stability": Only for one-level standards. A PDF visualizing the stability of the analytes in reference standards. Accepted deviations are displayed as dashed lines. 
- "IS_metric_plots": A scatter plot visualizing the stability of the internal standards over all samples. Accepted deviations are displayed as dashed lines. 
