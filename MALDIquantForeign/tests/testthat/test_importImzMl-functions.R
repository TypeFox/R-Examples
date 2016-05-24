context("importImzMl")

test_that("importImzMl continuous", {
  expect_error(MALDIquantForeign:::.importImzMl("tmp.tmp"))

  path <- normalizePath(system.file(
    file.path("exampledata", "tiny_continuous.imzML"),
    package="MALDIquantForeign"))

  s <- MALDIquantForeign:::.importImzMl(path, verbose=FALSE)

  expect_equal(s, import(path, verbose=FALSE))
  expect_equal(s, importImzMl(path, verbose=FALSE))
  expect_equal(s, import(path, type="imzML", verbose=FALSE))

  expect_equal(mass(s[[1]]), 1:5)
  expect_equal(intensity(s[[1]]), 6:10)
  expect_equal(basename(metaData(s[[1]])$file), "tiny_continuous.imzML")

  expect_equal(mass(s[[2]]), 1:5)
  expect_equal(intensity(s[[2]]), 10:6)
  expect_equal(basename(metaData(s[[2]])$file), "tiny_continuous.imzML")
})

test_that("importImzMl processed", {
  path <- normalizePath(system.file(
    file.path("exampledata", "tiny_processed.imzML"),
    package="MALDIquantForeign"))

  s <- MALDIquantForeign:::.importImzMl(path, verbose=FALSE)

  expect_equal(s, import(path, verbose=FALSE))
  expect_equal(s, importImzMl(path, verbose=FALSE))
  expect_equal(s, import(path, type="imzML", verbose=FALSE))

  expect_equal(mass(s[[1]]), 1:5)
  expect_equal(intensity(s[[1]]), 6:10)
  expect_equal(basename(metaData(s[[1]])$file), "tiny_processed.imzML")

  expect_equal(mass(s[[2]]), 6:10)
  expect_equal(intensity(s[[2]]), 10:6)
  expect_equal(basename(metaData(s[[2]])$file), "tiny_processed.imzML")
})

test_that("importImzMl coordinates", {
  path <- normalizePath(system.file(
    file.path("exampledata", "tiny_continuous.imzML"),
    package="MALDIquantForeign"))

  s <- MALDIquantForeign:::.importImzMl(path, verbose=FALSE)

  expect_equal(s, importImzMl(path, coordinates=cbind(1:2, c(1, 1)),
                              verbose=FALSE))
  expect_equal(s[1], importImzMl(path, coordinates=cbind(1, 1),
                                 verbose=FALSE))
  expect_equal(s[[2]], importImzMl(path, coordinates=cbind(2, 1),
                                   verbose=FALSE)[[1]])

  expect_error(importImzMl(path, coordinates=3, verbose=FALSE),
               "The .*coordinates.* argument has to be a matrix with two columns")
  expect_warning(importImzMl(path, coordinates=cbind(1:3, 2:0), verbose=FALSE),
                 "The following rows contain invalid coordinates: 1, 3")
})

