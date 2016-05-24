context("importMsd")

test_that("importMsd", {
  expect_error(MALDIquantForeign:::.importMsd("tmp.tmp"))

  path <- normalizePath(system.file(
    file.path("exampledata", "tiny1.msd"),
    package="MALDIquantForeign"))

  s <- MALDIquantForeign:::.importMsd(path, verbose=FALSE)

  expect_equal(s, import(path, verbose=FALSE))
  expect_equal(s, importMsd(path, verbose=FALSE))
  expect_equal(s, import(path, type="msd", verbose=FALSE))

  expect_true(isMassSpectrum(s[[1]]))
  expect_equal(mass(s[[1]]), 1:5)
  expect_equal(intensity(s[[1]]), 6:10)
  expect_equal(basename(metaData(s[[1]])$file), "tiny1.msd")
})
