library(testthat)

context("Validate")


test_that("validate_filename_output -good", {   
  path <- base::file.path(devtools::inst(name="IalsaSynthesis"), "test_data/2015-portland")
  
  invisible_1 <- validate_filename_output("b1_female_ae_muscle_fluid_grip_trailsb.out", path)
  invisible_2 <- validate_filename_output("b1_male_aehplus_muscle_memory_grip_logicalmemory.out", path)
  invisible_3 <- validate_filename_output("u1_male_aehplus_muscle_noCog_hand_noCogSpec.out", path)
  
  expect_true(invisible_1, "The first file name should be validated correctly.")
  expect_true(invisible_2, "The second file name should be validated correctly.")
  expect_true(invisible_3, "The third file name should be validated correctly.")
})

test_that("validate_filename_output -missing", {   
  path <- base::file.path(devtools::inst(name="IalsaSynthesis"), "test_data/2015-portland")
  expected_error_regex <- "The output file was not found at .+"
  expect_error(
    regexp = expected_error_regex,
    validate_filename_output("no_one_is_here.out", path)
  )
})

test_that("validate_filename_output -bad extension", {   
  path <- base::file.path(devtools::inst(name="IalsaSynthesis"), "test_data/2015-portland")
  expected_error_regex <- "The output file extension `bad` did not match the expected value of `out`."
  expect_error(
    regexp = expected_error_regex,
    validate_filename_output("extension.bad", path)
  )
})
