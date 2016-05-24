
context("mark_my_assignment")

test_that(desc="mark_my_assignment()",{
  suppressMessages(set_assignment(file.path(system.file(package = "markmyassignment"), "extdata/example_assignment01.yml")))
  
  source(file.path(system.file(package = "markmyassignment"), "extdata/example_lab_file.R"))
  
  expect_is(capture.output(mark_my_assignment()), "character")
  expect_error(capture.output(mark_my_assignment(reporter = "non_existing_reporter")))
  expect_is(mark_my_assignment(quiet = TRUE), "testthat_results")
  expect_is(mark_my_assignment(tasks = "task1", quiet = TRUE), "testthat_results")
  expect_equal(length(mark_my_assignment(tasks = "task1", quiet = TRUE)), 2)
  expect_is(mark_my_assignment(tasks = c("task1", "task2"), quiet = TRUE), "testthat_results")
  expect_is(mark_my_assignment(force_get_tests = TRUE, quiet = TRUE), "testthat_results")
  expect_warning(mark_my_assignment(mark_file = file.path(system.file(package = "markmyassignment"), "extdata/example_lab_file.R"), quiet = TRUE))
  
})


test_that(desc="Assertions on arguments in mark_my_assignment()",{
  
  suppressMessages(set_assignment(file.path(system.file(package = "markmyassignment"), "extdata/example_assignment01.yml")))
  source(file.path(system.file(package = "markmyassignment"), "extdata/example_lab_file.R"))
  
  expect_warning(mark_my_assignment(tasks = "no such task", quiet = TRUE))
  expect_warning(mark_my_assignment(tasks = c("task1", "no such task", "task2"), quiet = TRUE))
  expect_error(mark_my_assignment(tasks = task2, quiet = TRUE))
  expect_error(mark_my_assignment(quiet = "TRUE"))
  expect_error(mark_my_assignment(force_get_tests = "TRUE"))
  
})


test_that(desc="mark_my_dir()",{
  test_assgn_file <- file.path(system.file(package = "markmyassignment"), "extdata/example_assignment01.yml")
  test_dir <- file.path(system.file(package = "markmyassignment"), "extdata/example_dir")
  capture_output(res_mark <- mark_my_dir(directory = test_dir, lab_file = test_assgn_file))
  expect_is(res_mark, class = "list")
  expect_equal(length(res_mark), 2)
  expect_is(res_mark[[1]], class = "testthat_results")
})
