context("importBrukerFlex")

test_that("importBrukerFlex", {
  expect_error(MALDIquantForeign:::.importBrukerFlex("tmp.tmp"))

  path <- system.file(
    file.path("exampledata", "brukerflex", "0_A1", "1", "1SLin", "fid"),
    package="MALDIquantForeign")
  s <- MALDIquantForeign:::.importBrukerFlex(path, verbose=FALSE)

  expect_equal(s, import(path, verbose=FALSE))
  expect_equal(s, importBrukerFlex(path, verbose=FALSE))
  expect_equal(s, import(path, type="fid", verbose=FALSE))

  expect_equal(trunc(mass(s[[1]])), 226:230)
  expect_equal(intensity(s[[1]]), 1:5)
  expect_equal(basename(metaData(s[[1]])$file), "fid")
  expect_equal(metaData(s[[1]])$laserShots, 100)
  expect_equal(metaData(s[[1]])$comments, paste0("TESTSAMPLE", 1:4))
})
