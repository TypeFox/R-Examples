library(testthat)

context("Extract")

path_1 <- base::file.path(devtools::inst(name="IalsaSynthesis"), "test_data/2015-portland/b1_female_ae_muscle_fluid_grip_trailsb.out")
path_2 <- base::file.path(devtools::inst(name="IalsaSynthesis"), "test_data/2015-portland/b1_male_aehplus_muscle_memory_grip_logicalmemory.out")
path_3 <- base::file.path(devtools::inst(name="IalsaSynthesis"), "test_data/2015-portland/u1_male_aehplus_muscle_noCog_hand_noCogSpec.out")
output_1 <- readr::read_file(path_1)
output_2 <- readr::read_file(path_2)
output_3 <- readr::read_file(path_3)


snippet_fit_1 <- "

MODEL FIT INFORMATION

Number of Free Parameters                       25

Loglikelihood

H0 Value                      -11183.730

Information Criteria

Akaike (AIC)                   22417.460
Bayesian (BIC)                 22526.536
Sample-Size Adjusted BIC       22447.171
(n* = (n + 2) / 24)



MODEL RESULTS

Two-Tailed"

test_that("Path Data File", {   
  expected_1 <- "C:\\Users\\Dragon Hat\\Desktop\\CbaMaster.csv"
  expected_2 <- "C:\\Users\\Andrea Zammit\\Desktop\\EASMaster.csv"
  expected_3 <- "\"C:\\Users\\wuche\\Dropbox\\IALSA\\Data\\HABC-9999.dta.dat\"" #Notice this one has a fishy extension and enclosing quotes.
  
  snippet_1 <- "        TITLE: m5, b1,trails, peak flow, LGM,ae Conditional,female
  
  DATA:  File = C:\\Users\\Dragon Hat\\Desktop\\CbaMaster.csv;
  
  
  VARIABLE: Names are
  ! demographics
  SubjectID	Sex	Ethnic	Caus	DemEver	Bagesq	deathage	DEMO23	Educyrs	Status	"
  
  observed_from_snippet_1 <- IalsaSynthesis::extract_output_filename(snippet_1)
  observed_from_file_1 <- IalsaSynthesis::extract_output_filename(output_1)
  observed_from_file_2 <- IalsaSynthesis::extract_output_filename(output_2)
  observed_from_file_3 <- IalsaSynthesis::extract_output_filename(output_3)
  
  expect_equal(observed_from_snippet_1, expected_1, "The data file path extracted from the snippet should be correct.")
  expect_equal(observed_from_file_1, expected_1, "The data file path extracted from the first output file should be correct.")
  expect_equal(observed_from_file_2, expected_2, "The data file path extracted from the second output file should be correct.")
  expect_equal(observed_from_file_3, expected_3, "The data file path extracted from the third output file should be correct.")
})

test_that("Free Parameter Count", {   
  tolerance <- 0
  expected_1 <- 25
  expected_2 <- 41
  expected_3 <- 18
  
  observed_from_snippet_1 <- IalsaSynthesis::extract_free_parameter_count(snippet_fit_1)
  observed_from_file_1 <- IalsaSynthesis::extract_free_parameter_count(output_1)
  observed_from_file_2 <- IalsaSynthesis::extract_free_parameter_count(output_2)
  observed_from_file_3 <- IalsaSynthesis::extract_free_parameter_count(output_3)
  
  expect_equal(observed_from_snippet_1, expected_1, info="The free parameter count extracted from the snippet should be correct.", tolerance=tolerance)
  expect_equal(observed_from_file_1, expected_1, info="The free parameter count extracted from the first output file should be correct.", tolerance=tolerance)
  expect_equal(observed_from_file_2, expected_2, info="The free parameter count extracted from the second output file should be correct.", tolerance=tolerance)
  expect_equal(observed_from_file_3, expected_3, info="The free parameter count extracted from the third output file should be correct.", tolerance=tolerance)
})

test_that("Loglikelihood", {   
  tolerance <- 0.001
  expected_1 <- -11183.730
  expected_2 <- -1767.351
  expected_3 <- -20217.511
  
  observed_from_snippet_1 <- IalsaSynthesis::extract_loglikelihood(snippet_fit_1)
  observed_from_file_1 <- IalsaSynthesis::extract_loglikelihood(output_1)
  observed_from_file_2 <- IalsaSynthesis::extract_loglikelihood(output_2)
  observed_from_file_3 <- IalsaSynthesis::extract_loglikelihood(output_3)
  
  expect_equal(observed_from_snippet_1, expected_1, info="The free parameter count extracted from the snippet should be correct.", tolerance=tolerance)
  expect_equal(observed_from_file_1, expected_1, info="The loglikelihood extracted from the first output file should be correct.", tolerance=tolerance)
  expect_equal(observed_from_file_2, expected_2, info="The loglikelihood extracted from the second output file should be correct.", tolerance=tolerance)
  expect_equal(observed_from_file_3, expected_3, info="The loglikelihood extracted from the third output file should be correct.", tolerance=tolerance)
})

test_that("Scaling Correction Factor", {   
  tolerance <- 0.001
  # expected_1 <- NA_real_
  # expected_2 <- NA_real_
  expected_3 <- 1.3467
  
  observed_from_snippet_1 <- IalsaSynthesis::extract_scaling_correction(snippet_fit_1)
  observed_from_file_1 <- IalsaSynthesis::extract_scaling_correction(output_1)
  observed_from_file_2 <- IalsaSynthesis::extract_scaling_correction(output_2)
  observed_from_file_3 <- IalsaSynthesis::extract_scaling_correction(output_3)
  
  expect_true(is.na(observed_from_snippet_1), info="The Scaling Correction Factor extracted from the snippet should be NA")
  expect_true(is.na(observed_from_file_1), info="The Scaling Correction Factor extracted from the first output file should be NA")
  expect_true(is.na(observed_from_file_2), info="The Scaling Correction Factor extracted from the second output file should be NA")
  expect_equal(observed_from_file_3, expected_3, info="The Scaling Correction Factor extracted from the third output file should be correct.", tolerance=tolerance)
})

test_that("AIC", {   
  tolerance <- 0.001
  expected_1 <- 22417.460
  expected_2 <- 3616.701
  expected_3 <- 40471.022
  
  observed_from_snippet_1 <- IalsaSynthesis::extract_aic(snippet_fit_1)
  observed_from_file_1 <- IalsaSynthesis::extract_aic(output_1)
  observed_from_file_2 <- IalsaSynthesis::extract_aic(output_2)
  observed_from_file_3 <- IalsaSynthesis::extract_aic(output_3)
  
  expect_equal(observed_from_snippet_1, expected_1, info="The AIC extracted from the snippet should be correct.", tolerance=tolerance)
  expect_equal(observed_from_file_1, expected_1, info="The AIC extracted from the first output file should be correct.", tolerance=tolerance)
  expect_equal(observed_from_file_2, expected_2, info="The AIC extracted from the second output file should be correct.", tolerance=tolerance)
  expect_equal(observed_from_file_3, expected_3, info="The AIC extracted from the third output file should be correct.", tolerance=tolerance)
})

test_that("BIC", {   
  tolerance <- 0.001
  expected_1 <- 22526.536
  expected_2 <- 3710.044
  expected_3 <- 40566.198
  
  observed_from_snippet_1 <- IalsaSynthesis::extract_bic(snippet_fit_1)
  observed_from_file_1 <- IalsaSynthesis::extract_bic(output_1)
  observed_from_file_2 <- IalsaSynthesis::extract_bic(output_2)
  observed_from_file_3 <- IalsaSynthesis::extract_bic(output_3)
  
  expect_equal(observed_from_snippet_1, expected_1, info="The BIC extracted from the snippet should be correct.", tolerance=tolerance)
  expect_equal(observed_from_file_1, expected_1, info="The BIC extracted from the first output file should be correct.", tolerance=tolerance)
  expect_equal(observed_from_file_2, expected_2, info="The BIC extracted from the second output file should be correct.", tolerance=tolerance)
  expect_equal(observed_from_file_3, expected_3, info="The BIC extracted from the third output file should be correct.", tolerance=tolerance)
})

test_that("BIC: Sample-Size Adjusted", {   
  tolerance <- 0.001
  expected_1 <- 22447.171
  expected_2 <- 3580.867
  expected_3 <- 40509.018
  
  observed_from_snippet_1 <- IalsaSynthesis::extract_bic_adjusted(snippet_fit_1)
  observed_from_file_1 <- IalsaSynthesis::extract_bic_adjusted(output_1)
  observed_from_file_2 <- IalsaSynthesis::extract_bic_adjusted(output_2)
  observed_from_file_3 <- IalsaSynthesis::extract_bic_adjusted(output_3)

  expect_equal(observed_from_snippet_1, expected_1, info="The Sample-Size Adjusted BIC extracted from the snippet should be correct.", tolerance=tolerance)
  expect_equal(observed_from_file_1, expected_1, info="The Sample-Size Adjusted BIC extracted from the first output file should be correct.", tolerance=tolerance)
  expect_equal(observed_from_file_2, expected_2, info="The Sample-Size Adjusted BIC extracted from the second output file should be correct.", tolerance=tolerance)
  expect_equal(observed_from_file_3, expected_3, info="The Sample-Size Adjusted BIC extracted from the third output file should be correct.", tolerance=tolerance)
})
