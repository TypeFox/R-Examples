# Test files
good_file <- "data/Blending_001L00XS4.txt"

context("Printing and previewing data")

test_that("preview levels", {
  eprime_log <- read_eprime(good_file)
  chunked <- FrameList(eprime_log)
  output <- capture.output(preview_levels(chunked))

  # This file has 5 levels
  expect_equal(length(output), 5 + 2)
  expect_equal(output[1], "Level Counts: ")
})

test_that("preview frames", {
  eprime_log <- read_eprime(good_file)
  chunked <- FrameList(eprime_log)
  output <- capture.output(preview_frames(chunked))

  # First chunk is the Header
  expect_true(output[2] == " Eprime.Level Running Procedure")
  expect_true(output[3] == "            1  Header    Header")
})

test_that("preview eprime combines level and frame previews", {
  eprime_log <- read_eprime(good_file)
  chunked <- FrameList(eprime_log)
  output <- capture.output(preview_eprime(chunked))
  output_levels <- capture.output(preview_levels(chunked))
  output_frames <- capture.output(preview_frames(chunked))

  expect_true(all(output == c(output_levels, output_frames)))
})


test_that("print an EprimeFrame", {
  eprime_log <- read_eprime(good_file)
  chunked <- FrameList(eprime_log)

  # Print the frame
  frame <- chunked[[1]]
  printed_frame <- capture.output(frame)

  # Check specific lines in the print function
  proc_line <- ' $ Procedure          : chr "Header"'
  classes <- ' - attr(*, "class")= chr [1:2] "EprimeFrame" "list"'
  expect_equal(printed_frame[5 + 1], proc_line)
  expect_equal(printed_frame[22 + 2], classes)
})

test_that("print a FrameList", {
  eprime_log <- read_eprime(good_file)
  chunked <- FrameList(eprime_log)

  # Print the frame
  printed_list <- capture.output(chunked)
  nlines <- length(printed_list)

  # Check specific lines in the print function
  l1 <- "List of 25"
  l2 <- " $ :List of 22"
  last_1 <- '  ..- attr(*, "class")= chr [1:2] "EprimeFrame" "list"'
  last_2 <- ' - attr(*, "class")= chr [1:2] "list" "FrameList"'

  expect_equal(printed_list[1], l1)
  expect_equal(printed_list[2], l2)
  expect_equal(printed_list[nlines - 1], last_1)
  expect_equal(printed_list[nlines], last_2)
})

