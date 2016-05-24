## Test functions that create lists from character vectors of colon-separated
## values.

# expect_identical lists to have the same order
expect_nearly_identical <- function(x, y, ...) {
  expect_identical(sort_names(x), sort_names(y), ...)
}
sort_names <- function(x) x[sort(names(x))]


context("EprimeFrame")

test_that("Default values", {
  # Empty EprimeFrame
  default_list <- structure(list(
    Eprime.LevelName = NA,
    Eprime.Level = NA,
    Eprime.Basename = NA,
    Eprime.FrameNumber = NA,
    Procedure = NA,
    Running = NA), class = c("EprimeFrame", "list"))

  expect_nearly_identical(as.EprimeFrame(list()), default_list)

  # It's not possible to have an NA Eprime.Level when constructing from a
  # character vector
  default_character <- default_list
  default_character$Eprime.Level <- 1
  expect_nearly_identical(EprimeFrame(character()), default_character)
})

test_that("Simple construction", {
  test_values <- c(
    "key: value",
    "question: answer",
    "garbage text")

  expected_list <- structure(list(
    Eprime.LevelName = NA,
    Eprime.Level = 1,
    Eprime.Basename = NA,
    Eprime.FrameNumber = NA,
    Procedure = NA,
    Running = NA,
    key = "value",
    question = "answer"), class = c("EprimeFrame", "list"))

  expect_nearly_identical(EprimeFrame(test_values), expected_list)
})

test_that("Running-related values are normalized", {
  keys_values <- c(
    "Running: Demo",
    "Demo: ExampleCode",
    "Demo.Cycle: 1",
    "Demo.Sample: 1",
    "Key: Value")

  normalized_running <- structure(list(
    Eprime.Level = 1,
    Eprime.LevelName = "Demo_ExampleCode",
    Eprime.Basename = NA,
    Eprime.FrameNumber = NA,
    Procedure = NA,
    Running = "Demo",
    Cycle = "1",
    Sample = "1",
    Key = "Value"), class = c("EprimeFrame", "list"))

  expect_nearly_identical(EprimeFrame(keys_values), normalized_running)
})


test_that("listify handles typical data", {
  example_lines <- c(
    "\t*** LogFrame Start ***",
    "\tProcedure: PracticeProc",
    "\tStimulus1: toothbrush",
    "\t*** LogFrame End ***")

  expected_example <- list(
    Procedure = "PracticeProc",
    Stimulus1 = "toothbrush")

  expect_identical(listify(example_lines), expected_example)
})

test_that("listify handles garbage values", {
  garbage_lines <- c("", NA,
    "\t*** LogFrame Start ***",
    "*** LogFrame End ***")

  usable_lines <- c(
    "Item: 1",
    "\tWhiteSpaceItem: Stuff  ",
    "AnotherItem: Two Colons: Yes",
    "NoValue: ",
    "Last Item: Whatever")

  expected <- list(
    Item = "1",
    WhiteSpaceItem = "Stuff",
    AnotherItem = "Two Colons: Yes",
    NoValue = "",
    `Last Item` = "Whatever")
  expect_identical(listify(c(garbage_lines, usable_lines)), expected)

})
