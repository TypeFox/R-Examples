
context("set_assignment")

test_that(desc="set_assignment()",{
  correct_url1 <- "https://raw.githubusercontent.com/MansMeg/markmyassignment/master/inst/extdata/example_assignment01.yml"
  correct_url2 <- "https://raw.githubusercontent.com/MansMeg/markmyassignment/master/inst/extdata/example_assignment02.yml"
  correct_url3 <- "https://raw.githubusercontent.com/MansMeg/markmyassignment/master/inst/extdata/example_assignment03.yml"
  correct_url4 <- "https://raw.githubusercontent.com/MansMeg/markmyassignment/master/inst/extdata/example_assignment04.yml"
  correct_url5 <- "https://raw.githubusercontent.com/MansMeg/markmyassignment/master/inst/extdata/example_assignment05.yml"
  wrong_url1 <- "https://raw.githubusercontent.com/MansMeg/markmyassignment/master/inst/extdata/example_lab_file.R"
  wrong_url2 <- "https://raw.githubusercontent.com/MansMeg/markmyassignment/master/inst/extdata/file_that_do_not_exist.R"
  super_wrong_path <- "XXX"
  correct_local1 <- file.path(system.file(package = "markmyassignment"), "extdata/example_assignment01.yml")
  correct_local2 <- file.path(system.file(package = "markmyassignment"), "extdata/example_assignment02.yml")
  wrong_local1 <- file.path(system.file(package = "markmyassignment"), "extdata/example_lab_file.R")
  wrong_local2 <- file.path(system.file(package = "markmyassignment"), "file_that_do_not_exist.R")
  
  expect_is(suppressMessages(set_assignment(correct_url1)), "character")
  expect_is(suppressMessages(set_assignment(path = correct_url2)), "character")
  expect_is(suppressMessages(set_assignment(path = correct_url4)), "character")
  expect_is(suppressMessages(set_assignment(correct_local1)), "character")
  expect_is(suppressMessages(set_assignment(correct_local2)), "character")
  expect_error(set_assignment(path = correct_url3))
  expect_error(set_assignment(path = correct_url5))
  expect_error(set_assignment(path = wrong_url1))
  expect_error(set_assignment(wrong_url2))
  expect_error(set_assignment(wrong_local1))
  expect_error(set_assignment(wrong_local2))
  expect_error(set_assignment(super_wrong_path))
})

test_that(desc="show_tasks()",{
  correct_local1 <- file.path(system.file(package = "markmyassignment"), "extdata/example_assignment01.yml")
  correct_local2 <- file.path(system.file(package = "markmyassignment"), "extdata/example_assignment02.yml")
  suppressMessages(set_assignment(correct_local1))
  expect_equal(show_tasks(), c("task1","task2"))
  suppressMessages(set_assignment(correct_local2))
  expect_equal(show_tasks(), c("task1","task2"))
})

test_that(desc="check_installed_packages()",{
  assgn_path <- file.path(system.file(package = "markmyassignment"), "extdata/example_assignment08_bad_pkgs.yml")
  expect_warning(set_assignment(path = assgn_path), regexp = "The following packages need to be installed and then loaded")
  assgn_path <- file.path(system.file(package = "markmyassignment"), "extdata/example_assignment07_pkgs.yml")
  expect_warning(set_assignment(path = assgn_path), regexp = "The following packages need to be loaded")
  library(codetools)
  expect_is(suppressMessages(set_assignment(assgn_path)), "character")
  detach(name = "package:codetools")
})
