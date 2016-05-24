
context("mark_my_file")

test_that(desc="mark_my_file()",{
  source_file <- file.path(system.file(package = "markmyassignment"), "extdata/example_lab_file.R")
  assignment_file <- file.path(system.file(package = "markmyassignment"), "extdata/example_assignment01.yml")
  
  expect_is(mark_my_file(mark_file = source_file, assignment_path = assignment_file, quiet = TRUE), "testthat_results")
  expect_is(mark_my_file(mark_file = file.path(system.file(package = "markmyassignment"), "extdata/example_lab_file.R"), assignment_path = assignment_file, quiet = TRUE), "testthat_results")  
  expect_is(mark_my_file(mark_file = file.path(system.file(package = "markmyassignment"), "extdata/example_lab_file_circular.R"), assignment_path = assignment_file, quiet = TRUE), "testthat_results")
  capture_output(expect_is(mark_my_file(mark_file = file.path(system.file(package = "markmyassignment"), "extdata/example_lab_file_messy.R"), assignment_path = assignment_file, quiet = TRUE), "testthat_results"))
  expect_is(capture.output(mark_my_file(mark_file = source_file, assignment_path = assignment_file)), "character")
  expect_is(mark_my_file(tasks = "task1", mark_file = source_file, assignment_path = assignment_file, quiet = TRUE), "testthat_results")
  expect_equal(length(mark_my_file(tasks = "task1", mark_file = source_file, assignment_path = assignment_file, quiet = TRUE)), 2)
  expect_is(mark_my_file(tasks = c("task1", "task2"), mark_file = source_file, assignment_path = assignment_file, quiet = TRUE), "testthat_results")
  expect_is(mark_my_file(mark_file = source_file, assignment_path = assignment_file, force_get_tests = TRUE, quiet = TRUE), "testthat_results")
})

test_that(desc="Assertions on arguments in mark_my_file()",{
  source_file <- file.path(system.file(package = "markmyassignment"), "extdata/example_lab_file.R")
  assignment_file <- file.path(system.file(package = "markmyassignment"), "extdata/example_assignment01.yml")
  
  expect_warning(mark_my_file(tasks = "no such task", mark_file = source_file, assignment_path = assignment_file, quiet = TRUE))
  expect_error(mark_my_file(mark_file = source_file, assignment_path = assignment_file, tasks = task2, quiet = TRUE))
  expect_error(mark_my_file(mark_file = source_file, assignment_path = "~/no such directory/no such file.yml", quiet = TRUE))
  expect_error(mark_my_file(quiet = "TRUE", mark_file = source_file, assignment_path = assignment_file))
  expect_error(mark_my_file(force_get_tests = "TRUE", mark_file = source_file, assignment_path = assignment_file))
  expect_error(mark_my_file(reporter = StudentReporter, mark_file = source_file, assignment_path = assignment_file))
  expect_error(mark_my_file(mark_file = "~/no such directory/no such file.R"))
})


test_that(desc="Load packages before running mark_my_file()",{
  skip("These tests need to be fixed")
  source_file <- file.path(system.file(package = "markmyassignment"), "extdata/example_lab_file.R")
  assignment_file <- file.path(system.file(package = "markmyassignment"), "extdata/example_assignment08_bad_pkgs.yml")
  expect_error(mark_my_file(mark_file = source_file, assignment_path = assignment_file), regexp = "The following packages need to be installed and then loaded")
  
  assignment_file <- file.path(system.file(package = "markmyassignment"), "extdata/example_assignment07_pkgs.yml")
  expect_error(mark_my_file(mark_file = source_file, assignment_path = assignment_file), regexp = "The following packages need to be loaded")
  library(codetools)
  expect_is(mark_my_file(mark_file = source_file, assignment_path = assignment_file, quiet = TRUE), "testthat_results")
  detach(name = "package:codetools")
  }
)
