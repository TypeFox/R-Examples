################################################################################
# test_infrastructure.R
# (c)2012 Andrew Redd 
# This is file part of the lint R package, a code style check package for R.
# 
# DESCRIPTION
# -----------
# This file contains the tests for checking the infrastructure and utilities.
# 
################################################################################
context("Infrastructure")
library(stringr)

test_that('find_region returns empty.find for empty.regions',{
    expect_that(find_region(character(0),,,), is_identical_to(empty.find))
})
test_that("with_default",{
  expect_that(with_default(NULL, T), is_true())
  expect_that(with_default(NULL, F), is_false())
  expect_that(with_default(NA, T), is_true())
  expect_that(with_default(NA, F), is_false())
  expect_that(with_default(NULL, "test"), is_identical_to("test"))
  expect_that(with_default(character(0), "test"), is_identical_to(character(0)))
  expect_that(with_default(T, F), is_true())
  expect_that(with_default(F, T), is_false())
  expect_that(with_default("test", T), is_identical_to("test"))
})
test_that("check_pattern", {
  lines = c(
    "123",
    "abc",
    "xyz")
  expect_identical(check_pattern(lines, "no match"), empty.find
    , 'no match gives and empty find result.')
  expect_that(check_pattern(lines, "123")$line1, is_equivalent_to(1L)
    , 'find the right line')
  expect_that(valid_find(check_pattern(lines, "123")), is_true()
    , 'return type is a valid find formatted data.frame')
})
test_that("dispatch_test", {
  file <- 
  check.file <- find_example("checks.R", package="lint")
  pd <-
  parse.data <- getParseData(parse(check.file, keep.source=TRUE))
  expect_that(
    dispatch_test(list(exclude.region=character(0)), check.file)
  , throws_error("Ill-formatted check."))

  lines = c(
    "123",
    "abc",
    "xyz")
  parse.data <- 
  pd <- getParseData(parse(text=paste(lines,'\n', collapse=''), keep.source=TRUE))
  test <- list(pattern='abc')
  expect_that(
      dispatch_test(test, , pd, lines, warning=T)
    , gives_warning('Lint: abc: found on lines 2'))
  expect_that(
      dispatch_test(list(pattern=perl("\\w{3}")), 
                   , parse.data=pd, lines=lines,  warning=T)
    , gives_warning('Lint: .*: found on lines 1, 2, 3'))
  expect_that(
      dispatch_test(list(pattern=perl("\\d{3}")), 
                   , parse.data=pd, lines=lines,  warning=T)
    , gives_warning('Lint: .*: found on lines 1'))
})
test_that("lint", {
    file <- 
    check.file <- find_example("checks.R", package="lint")
        
    lint.tests <- list(
          spacing.twobeforecomments  = spacing.twobeforecomments
        , spacing.spacearoundequals  = spacing.spacearoundequals
        # , spacing.indentation.notabs = spacing.indentation.notabs
        , spacing.linelength.80      = spacing.linelength.80)
    lint(check.file, lint.tests)
})

{# span_intersect
span.1119 <- make_ex_span(1, 1, 1, 9)
span.1114 <- make_ex_span(1, 1, 1, 4)
span.1517 <- make_ex_span(1, 5, 1, 7)
span.2129 <- make_ex_span(2, 1, 2, 9)
test_that('span_intersect returns empty.find when given empty arguments', {
    expect_identical(span_intersect(span.1119, empty.find), empty.find)
    expect_identical(span_intersect(empty.find, span.1119), empty.find)
})
test_that('span_intersect can identify nested spans.', {
    expect_identical(span_intersect(span.1119, span.1517), span.1517)
    expect_identical(span_intersect(span.1517, span.1119), span.1517)
})
test_that('span_intersect returns empty when nonintersecting spans are given', {
    expect_identical(span_intersect(span.1114, span.1517), empty.find)
    expect_identical(span_intersect(span.1119, span.2129), empty.find)
})
test_that('span_intersect can recognize multiple overlaps', {
    spans.A <- rbind(span.1114, span.1517)
    expect_identical(span_intersect(spans.A, span.1119), spans.A)
    expect_identical(span_intersect(span.1119, spans.A), spans.A)
})
test_that('span_intersect can correctly do overlaps', {
    span.1117 <- make_ex_span(1, 1, 1, 7)
    span.1519 <- make_ex_span(1, 5, 1, 9)
    span.1527 <- make_ex_span(1, 5, 2, 7)
    span.2127 <- make_ex_span(2, 1, 2, 7)
    expect_identical(span_intersect(span.1117, span.1519), span.1517)
    expect_identical(span_intersect(span.1519, span.1117), span.1517)
    expect_identical(span_intersect(span.1527, span.2129), span.2127)
    expect_identical(span_intersect(span.2129, span.1527), span.2127)
})
}

test_that("format_problem_lines",{
    expect_identical(format_problem_lines(1:3, 5), '1, 2, 3')
    expect_identical(format_problem_lines(1:5, 5), '1, 2, 3, 4, 5')
    expect_identical(format_problem_lines(1:6, 5), '1, 2, 3, 4, 5, +1 more.')
})
