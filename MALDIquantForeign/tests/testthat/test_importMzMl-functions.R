context("importMzMl")

test_that("importMzMl", {
  expect_error(MALDIquantForeign:::.importMzMl("tmp.tmp"))

  path <- normalizePath(system.file(
    file.path("exampledata", "tiny1.mzML1.1.mzML"),
    package="MALDIquantForeign"))

  s <- MALDIquantForeign:::.importMzMl(path, verbose=FALSE)

  expect_equal(s, import(path, verbose=FALSE))
  expect_equal(s, importMzMl(path, verbose=FALSE))
  expect_equal(s, import(path, type="mzML", verbose=FALSE))

  expect_true(isMassSpectrum(s[[1]]))
  expect_equal(mass(s[[1]]), 1:5)
  expect_equal(intensity(s[[1]]), 6:10)
  expect_equal(basename(metaData(s[[1]])$file), "tiny1.mzML1.1.mzML")

  expect_true(isMassSpectrum(s[[2]]))
  expect_equal(mass(s[[2]]), 1:5)
  expect_equal(intensity(s[[2]]), 10:6)
  expect_equal(basename(metaData(s[[2]])$file), "tiny1.mzML1.1.mzML")
})

test_that("importMzMl compressed", {
  expect_error(MALDIquantForeign:::.importMzMl("tmp.tmp"))

  path <- normalizePath(system.file(
    file.path("exampledata", "tiny1-compressed.mzML1.1.mzML"),
    package="MALDIquantForeign"))

  s <- MALDIquantForeign:::.importMzMl(path, verbose=FALSE)

  expect_equal(s, import(path, verbose=FALSE))
  expect_equal(s, importMzMl(path, verbose=FALSE))
  expect_equal(s, import(path, type="mzML", verbose=FALSE))

  expect_true(isMassSpectrum(s[[1]]))
  expect_equal(mass(s[[1]]), 1:5)
  expect_equal(intensity(s[[1]]), 6:10)
  expect_equal(basename(metaData(s[[1]])$file), "tiny1-compressed.mzML1.1.mzML")

  expect_true(isMassSpectrum(s[[2]]))
  expect_equal(mass(s[[2]]), 1:5)
  expect_equal(intensity(s[[2]]), 10:6)
  expect_equal(basename(metaData(s[[2]])$file), "tiny1-compressed.mzML1.1.mzML")
})

test_that("importMzMl centroided", {
  expect_error(MALDIquantForeign:::.importMzMl("tmp.tmp"))

  path <- normalizePath(system.file(
    file.path("exampledata", "tiny1-centroided.mzML1.1.mzML"),
    package="MALDIquantForeign"))

  p <- MALDIquantForeign:::.importMzMl(path, centroided=TRUE, verbose=FALSE)

  expect_equal(p, import(path, centroided=TRUE, verbose=FALSE))
  expect_equal(p, importMzMl(path, centroided=TRUE, verbose=FALSE))
  expect_equal(p, import(path, type="mzML", centroided=TRUE, verbose=FALSE))

  expect_true(isMassPeaks(p[[1]]))
  expect_equal(mass(p[[1]]), 1:5)
  expect_equal(intensity(p[[1]]), 6:10)
  expect_equal(basename(metaData(p[[1]])$file), "tiny1-centroided.mzML1.1.mzML")

  expect_true(isMassPeaks(p[[2]]))
  expect_equal(mass(p[[2]]), 1:5)
  expect_equal(intensity(p[[2]]), 10:6)
  expect_equal(basename(metaData(p[[2]])$file), "tiny1-centroided.mzML1.1.mzML")

  ## overwrite default arguments
  path <- normalizePath(system.file(
    file.path("exampledata", "tiny1.mzML1.1.mzML"),
    package="MALDIquantForeign"))
  expect_true(all(sapply(MALDIquantForeign:::.importMzMl(path, centroided=TRUE),
                         isMassPeaks)))
})
