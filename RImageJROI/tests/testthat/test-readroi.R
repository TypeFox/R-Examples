path <- file.path(system.file(package = "RImageJROI"), "extdata", "ijroi")

context("Reading a rectangle")

test_that("dimensions are correct", {
  r <- read.ijroi(file.path(path, "rect.roi"))
  expect_that(r$left,   equals(31))
  expect_that(r$top,    equals(6))
  expect_that(r$right,  equals(94))
  expect_that(r$bottom, equals(25))
  expect_that(r$width,  equals(63))
  expect_that(r$height, equals(19))
})

context("Reading a polygon")

test_that("dimensions are correct", {
  coords <- structure(c(1, 31, 31, 23, 23, 31, 31, 1, 1, 9, 9, 1, 0, 0, 2, 
7, 18, 23, 25, 25, 23, 18, 7, 2), .Dim = c(12L, 2L))
  colnames(coords) <- c("x", "y")
  r <- read.ijroi(file.path(path, "polygon.roi"))
  expect_that(r$coords, equals(coords))
})

context("Reading an oval")

test_that("dimensions are correct", {
  r <- read.ijroi(file.path(path, "oval.roi"))
  expect_that(r$left,   equals(161))
  expect_that(r$top,    equals(6))
  expect_that(r$right,  equals(192))
  expect_that(r$bottom, equals(26))
  expect_that(r$width,  equals(31))
  expect_that(r$height, equals(20))
})

context("Reading test file with options")
test_that("options are read", {
  ## This file has OUTLINE and DOUBLE_HEADED set
  r <- read.ijroi(file.path(path, "arrow-outline-doubleheaded.roi"))
  expect_that(r$options, equals(6))
  expect_that(r$doubleHeaded, equals(TRUE))
  expect_that(r$outline, equals(TRUE))
})

