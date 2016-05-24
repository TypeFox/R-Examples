## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  eval = FALSE
  
)

## ---- eval = TRUE, echo = FALSE------------------------------------------
example_bin_counts <- structure(c(308, 243, 0, 624, 448, 0, 247, 163, 303, 263, 189, 
637, 249, 174, 266, 216, 270, 262, 287, 257, 274, 213, 234, 244, 
157, 269, 245, 240, 262, 142, 217, 209, 148, 225, 264, 229, 239, 
231, 142, 243, 99, 225, 251, 241, 205, 324, 88, 255, 146, 285, 
258, 251, 241, 355, 264, 259, 236, 259, 242, 216, 244, 383, 230, 
261), .Dim = c(8L, 8L), .Dimnames = list(c(" 1", "2", "3", "4", 
"5", "6", "7", "8"), c("1", "2", "3", "4", "5", "6", "7", "8"
)))

## ---- eval = TRUE, echo = FALSE------------------------------------------
example_regression_stats <- structure(c("-0.171248675025602", "0.00615122686158584", "Practical_CV", 
"0.583158688871931", "11  2  5  12", "1.00050469190991", "0.00335713985926308", 
"-0.799303637284865", "0.00599131304613039", "Practical_CV", 
"0.503865450973465", "10  19  17  16", "1.00011356353248", "0.00394053718934269"
), .Dim = c(7L, 2L), .Dimnames = list(c("Z_score_sample", "CV", 
"cv_types", "P_value_shapiro", "Predictor_chromosomes", "Mean_test_set", 
"CV_train_set"), c("Prediction_set_1", "Prediction_set_2")))

## ---- eval = TRUE, echo = FALSE------------------------------------------
example_bin_counts_sex <- structure(c(0, 0, 160, 0, 6, 0, 266, 0, 127, 0, 85, 0, 182, 0, 
235, 0), .Dim = c(2L, 8L), .Dimnames = list(c("X", "Y"), c("1", 
"2", "3", "4", "5", "6", "7", "8")))

## ---- eval = TRUE, echo = FALSE------------------------------------------
example_match_report <- structure(c(6.76557504762928e-07, 2.15862500903867e-06, 3.10100948006486e-06, 
3.29198583437393e-06, 3.58809416090935e-06, 3.83885934221371e-06, 
4.43239735020721e-06, 7.8916102881966e-06, 9.73718696931918e-06, 
1.12605046722957e-05), .Dim = c(10L, 1L), .Dimnames = list(c("Sample 6", 
"Sample 8", "Sample 1", "Sample 10", "Sample 9", "Sample 3", "Sample 5", 
"Sample 2", "Sample 7", "Sample 4"), "Sum_of_squares"))

## ---- eval=TRUE, echo = FALSE, results = "asis"--------------------------
pander::pandoc.table(example_bin_counts, caption = "Table 1: Example data bin counts. The rows represent chromosomes and the columns bins", include.rownames = T)

## ---- eval=TRUE, echo = FALSE, results = "asis"--------------------------
pander::pandoc.table(example_bin_counts_sex, caption = "Table 2: Example data bin counts sex chromosomes. The rows represent chromosomes and the columns bins", include.rownames = T)

## ------------------------------------------------------------------------
#  bin_bam_sample(bam_filepath, do_sort = F, separate_strands = F, custom_name = NULL)

## ------------------------------------------------------------------------
#  gc_correct(nipt_object, method, include_XY = F, span = 0.75)

## ------------------------------------------------------------------------
#  #Correct NIPTSample object using LOESS method
#  loess_corrected_sample <- gc_correct(nipt_object = sample_of_interest, method = "LOESS",
#                                       include_XY = F, span = 0.75)
#  #Correct NIPTControlGroup object using bin method
#  gc_bin_corrected_control_group <- gc_correct(nipt_object = control_group, method = "bin",
#                                               include_XY = T)
#  

## ------------------------------------------------------------------------
#  matchcontrolgroup(nipt_sample, nipt_control_group, mode, n_of_samples,
#    include_chromosomes = NULL, exclude_chromosomes = NULL)

## ---- eval=TRUE, echo = FALSE, results = "asis"--------------------------
pander::pandoc.table(example_match_report , caption = "Table 3: Example match_control_group mode 'subset'", include.rownames = T)

## ------------------------------------------------------------------------
#  #Run matchcontrolgroup in mode "report"
#  scores_control_group <- matchcontrolgroup(nipt_sample = gc_LOESS_corrected_sample,
#                                            nipt_control_group = gc_LOESS_corrected_control_group,
#                                            mode = "report", include_chromosomes = c(13,18))
#  #Run matchcontrolgroup in mode "subset" and select 50 best matching samples
#  subset_loess_corrected_control_group <- matchcontrolgroup(nipt_sample = gc_LOESS_corrected_sample,
#                                                            nipt_control_group = loess_corrected_control_group,
#                                                            mode = "subset", n_of_samples = 50)
#  

## ------------------------------------------------------------------------
#  chi_correct(nipt_sample, nipt_control_group, chi_cutoff = 3.5, include_XY = F)

## ------------------------------------------------------------------------
#  #Apply chi-squared based variation reduction method
#  chi_corrected_data <- chicorrect(nipt_sample = gc_LOESS_corrected_sample,
#                                   nipt_control_group = subset_loess_corrected_control_group)
#  #Extract sample and control group
#  loess_chi_corrected_sample <- chi_corrected_data$sample
#  subset_loess_chi_corrected_control_group <- chi_corrected_data$control_group

## ------------------------------------------------------------------------
#  calculate_z_score(nipt_sample, nipt_control_group, chromo_focus)

## ------------------------------------------------------------------------
#  #Calculate Z score for chromosome 13
#  z_score_result_13 <- calculate_z_score(nipt_sample = loess_chi_corrected_sample,
#                                         nipt_control_group = subset_loess_chi_corrected_control_group,
#                                         chromo_focus = 13)

## ------------------------------------------------------------------------
#  perform_regression(nipt_sample, nipt_control_group, chromo_focus,
#    n_models = 4, n_predictors = 4, exclude_chromosomes = NULL,
#    include_chromosomes = NULL, use_test_train_set = T,
#    size_of_train_set = 0.6, overdispersion_rate = 1.15,
#    force_practical_vc = F)

## ---- eval=TRUE, echo = FALSE, results = "asis"--------------------------
pander::pandoc.table(example_regression_stats, caption = "Table 2: Example regression results and statistics", include.rownames = T)

## ------------------------------------------------------------------------
#  prepare_ncv(nipt_control_group, chr_focus, max_elements,
#    exclude_chromosomes = NULL, include_chromosomes = NULL,
#    use_test_train_set = T, size_of_train_set = 0.6)

## ------------------------------------------------------------------------
#  calculate_ncv_score(nipt_sample, ncv_template)

## ------------------------------------------------------------------------
#  as_control_group(nipt_samples)

## ------------------------------------------------------------------------
#  add_samples_controlgroup(nipt_control_group, samples_to_add)

## ------------------------------------------------------------------------
#  remove_duplicates_controlgroup(nipt_control_group)

## ------------------------------------------------------------------------
#  remove_sample_controlgroup(samplename, nipt_control_group)

## ------------------------------------------------------------------------
#  diagnose_control_group(nipt_control_group)

## ------------------------------------------------------------------------
#  
#  list_NIPT_samples <- lapply(X = bam_filepaths, bin_bam_sample, do_sort = FALSE,
#                               separate_strands = T, custom_name = NULL)
#  

## ------------------------------------------------------------------------
#  control_group  <- as_control_group(nipt_samples = lapply(X = bam_filepaths, bin_bam_sample,
#                                                           do_sort = F, separate_strands = T))

## ------------------------------------------------------------------------
#  control_group  <- as_control_group(nipt_samples = mapply(FUN = bin_bam_sample, bam_filepaths,
#                                                           custom_name = names, SIMPLIFY = FALSE))

## ------------------------------------------------------------------------
#  library(NIPTeR)
#  #Retrieve filenames
#  bam_filepaths <- list.files(path = "/Path/to/bamfiles/", pattern = ".bam", full.names = T)
#  #Load files and convert to control group
#  control_group  <- as_control_group(nipt_samples = lapply(X = bam_filepaths, bin_bam_sample,
#                                                           do_sort = F, separate_strands = FALSE))
#  #Save control group for later
#  saveRDS(object = control_group, file = "/Path/to/directory/control_group.rds")

## ------------------------------------------------------------------------
#  library(NIPTeR)
#  #Gather all bam filepaths in a vector. Corresponds to 1a in figure
#  bam_filepaths <- list.files(path = "/Path/to/bamfiles/", pattern = ".bam", full.names = T)
#  #Load all bam files using lapply and feed the results to function as_control_group,
#  #converting the NIPTSamples to a NIPTControlGroup object. Corresponds to 1b in figure
#  control_group  <- as_control_group(nipt_samples = lapply(X = bam_filepaths, bin_bam_sample,
#                                                           do_sort = F, separate_strands = F))
#  #apply a gc LOESS correction to all samples. Since this can take up to 30 seconds
#  #sample, doing this once for a control group and store the results can be rewarding
#  #in terms of analysis time. Corresponds to 2a in figure
#  gc_control_group <- gc_correct(nipt_object = control_group, method = "LOESS")
#  #Retrieve control group diagnostics. Corresponds with 3a in figure
#  control_group_diagnostics <- diagnose_control_group(nipt_control_group = control_group)
#  #Retrieve samplenames with an abberant Z score for any chromosome and remove these samples
#  #from the control group. Corresponds with 3b in figure
#  abberant_sample_names <- unique(control_group_diagnostics$abberant_scores$Sample_name)
#  for (i in 1:length(abberant_sample_names)){
#    control_group <- remove_sample_controlgroup(samplename = abberant_sample_names[i],
#                                                nipt_control_group = control_group)
#  }
#  #Save the gc_corrected control groups to disk. Corresponds to 4a in figure
#  saveRDS(object = control_group, file = "/path/to/controlgroup.rds")

## ------------------------------------------------------------------------
#  library(NIPTeR)
#  #Load sample. Corresponds with 1a in figure
#  sample_of_interest <- bin_bam_sample(bam_filepath = "/Path/to/bamfile.bam", separate_strands = T)
#  #Load control group. Corresponds with 1b in figure
#  control_group <- readRDS("/Path/to/control_group.rds")
#  #Perform a chi-square based variation reduction on new trimmed control group and sample and
#  #extract data. Corresponds with 2a in figure
#  chi_data <- chi_correct(nipt_sample = sample_of_interest, nipt_control_group = control_group)
#  sample_of_interest <- chi_data$sample
#  control_group <- chi_data$control_group
#  #Perform regression for chromosome 21 with default settings, so:
#  #Create 4 models with 4 predictors each
#  #All chromosomes are potential predictors except the potential trisomic chromosomes 13, 18 and 21
#  #Use a test and train set where the size of the train set is 60% of the control group
#  #Assume at least 15% overdispersion in the data
#  #Dont force practical CV, so if the CV is below 1.15 * standard error of the mean use this as VC
#  #Corresponds with 3a in figure
#  regression_score_21 <- perform_regression(nipt_sample = sample_of_interest,
#                                            nipt_control_group = control_group, chromo_focus = 21)
#  
#  
#  ###       Gather relevant data from the objects on the workspace     #######

## ------------------------------------------------------------------------
#  library(NIPTeR)
#  #Load sample. Corresponds with 1a in figure
#  sample_of_interest <- bin_bam_sample(bam_filepath = "/Path/to/bamfile.bam")
#  #Load control group. Corresponds with 1b in figure
#  control_group <- readRDS("/Path/to/control_group.rds")
#  #Peform a GC correction type bin for the sample and the controlgroup
#  #Corresponds with 2a in figure
#  sample_of_interest <- gc_correct(nipt_object = sample_of_interest, method = "bin")
#  control_group <- gc_correct(nipt_object = control_group, method = "bin")
#  #Trim control group by only selecting 80% of best matching control samples
#  #Corresponds with 2b in figure
#  control_group <- match_control_group(nipt_sample = sample_of_interest, nipt_control_group = control_group,
#                                       mode = "subset", n_of_samples = round(length(control_group$samples) *.8,
#                                                                             digits = 0))
#  #Perform a chi-square based variation reduction on new trimmed control group and sample and
#  #extract data. Corresponds with 2c in figure
#  chi_data <- chi_correct(nipt_sample = sample_of_interest, nipt_control_group = control_group)
#  sample_of_interest <- chi_data$sample
#  control_group <- chi_data$control_group
#  #Retrieve control group diagnostics. Corresponds with 3a in figure
#  control_group_diagnostics <- diagnose_control_group(nipt_control_group = control_group)
#  #Retrieve samplenames with an abberant Z score for any chromosome and remove these samples
#  #from the control group. Corresponds with 3b in figure
#  abberant_sample_names <- unique(control_group_diagnostics$abberant_scores$Sample_name)
#  for (i in 1:length(abberant_sample_names)){
#    control_group <- remove_sample_controlgroup(samplename = abberant_sample_names[i],
#                                                nipt_control_group = control_group)
#  }
#  #Calculate Z score from chromosomes 13, 18 and 21. Corresponds with 4a in figure
#  z_score_13 <- calculate_z_score(nipt_sample = sample_of_interest,
#                                  nipt_control_group = control_group, chromo_focus = 13)
#  z_score_18 <- calculate_z_score(nipt_sample = sample_of_interest,
#                                  nipt_control_group = control_group, chromo_focus = 18)
#  z_score_21 <- calculate_z_score(nipt_sample = sample_of_interest,
#                                  nipt_control_group = control_group, chromo_focus = 21)
#  #Perform regression for all potential trisomic chromosomes with default settings, so:
#  #Create 4 models for every potential trisomic chromosome with 4 predictors each
#  #All chromosomes are potential predictors except the potential trisomic chromosomes 13, 18 and 21
#  #Use a test and train set where the size of the train set is 60% of the control group
#  #Assume at least 15% overdispersion in the data
#  #Dont force practical CV, so if the CV is below 1.15 * standard error of the mean use this as VC
#  #Corresponds with 4c in figure
#  regression_score_13 <- perform_regression(nipt_sample = sample_of_interest,
#                                            nipt_control_group = control_group, chromo_focus = 13)
#  regression_score_18 <- perform_regression(nipt_sample = sample_of_interest,
#                                            nipt_control_group = control_group, chromo_focus = 18)
#  regression_score_21 <- perform_regression(nipt_sample = sample_of_interest,
#                                            nipt_control_group = control_group, chromo_focus = 21)
#  #Get NCVTemplates for all potential trisomic chromosomes with max 9 denominators and default settings, so:
#  #All autosomals chromosomes are potential predictors, except the potential trisomic chromosomes 13, 18 and 21
#  #Use a test and train set where the size of the train set is 60% of the control group
#  #Corresponds with 4c in figure
#  new_ncv_template_13 <- prepare_ncv(nipt_control_group = control_group, chr_focus = 13, max_elements = 9)
#  new_ncv_template_18 <- prepare_ncv(nipt_control_group = control_group, chr_focus = 18, max_elements = 9)
#  new_ncv_template_21 <- prepare_ncv(nipt_control_group = control_group, chr_focus = 21, max_elements = 9)
#  #Use the NCVTemplates to get NCV scores for the sample of interest
#  #Corresponds with 4d in figure
#  ncv_score_13 <- calculate_ncv_score(nipt_sample = sample_of_interest, ncv_template = new_ncv_template_13)
#  ncv_score_18 <- calculate_ncv_score(nipt_sample = sample_of_interest, ncv_template = new_ncv_template_18)
#  ncv_score_21 <- calculate_ncv_score(nipt_sample = sample_of_interest, ncv_template = new_ncv_template_21)
#  
#  
#  ###       Gather relevant data from the objects on the workspace     #######
#  
#  
#  

