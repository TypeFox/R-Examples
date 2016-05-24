# Helpers
first <- function(...) head(..., n = 1)
last <- function(...) tail(..., n = 1)
first_line_header <- function(x) first(x) == "*** Header Start ***"
count_blanks <- function(xs) sum(stringr::str_detect(xs, "^$"))

# Test files
good_file <- "data/MINP_001L00XS1.txt"
bad_encoding <- "data/Blending_001L00XS4.txt"
crashed_experiment <- "data/MP_Block1_001P00XA1.txt"
no_footer <- "data/Coartic_Block1_001P00XS1.txt"
not_an_eprime_file <- "data/not_an_eprime_file.txt"


context("Reading standard files")

test_that("Load well-formed data", {
  eprime_log <- read_eprime(good_file, remove_clock = FALSE)
  expect_equal(first(eprime_log), "*** Header Start ***")
  expect_equal(eprime_log[29], "\t*** LogFrame Start ***")
  expect_equal(last(eprime_log), "*** LogFrame End ***")

  # Removing the clock removes 2 lines, one of which is in the header
  no_clock <- read_eprime(good_file, remove_clock = TRUE)
  expect_equal(no_clock[28], "\t*** LogFrame Start ***")
  expect_equal(length(no_clock) + 2, length(eprime_log))
})


context("Reading non-standard files")

test_that("load file with unexpected encoding", {
  # no warnings
  warnings <- evaluate_promise(read_eprime(bad_encoding))$warnings
  expect_equal(length(warnings), 0)

  # first line header and no blank lines
  bad_encoding_lines <- read_eprime(bad_encoding)
  expect_true(first_line_header(bad_encoding_lines))
  expect_less_than(count_blanks(bad_encoding_lines), 1)
})

test_that("load file from a crashed experiment", {
  # no warnings
  warnings <- evaluate_promise(read_eprime(crashed_experiment))$warnings
  expect_equal(length(warnings), 0)

  # first line header and no blank lines
  crashed_experiment_lines <- read_eprime(crashed_experiment)
  expect_true(first_line_header(crashed_experiment_lines))
  expect_less_than(count_blanks(crashed_experiment_lines), 1)
})

test_that("non-eprime file raises warning", {
  expect_warning(read_eprime(not_an_eprime_file), "not an Eprime txt file")
})

test_that("non-eprime file uses dummy text", {
  # first line header
  not_an_eprime_file_lines <- suppressWarnings(read_eprime(not_an_eprime_file))
  expect_true(first_line_header(not_an_eprime_file_lines))
})
