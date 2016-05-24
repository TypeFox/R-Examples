context("testChecksum")


test_that("testChecksum-functions", {
  f <- normalizePath(system.file(
    file.path("exampledata", "ascii.txt"),
    package="MALDIquantForeign"))
  md5 <- "9274dd34d675950326a222a952309a17"
  sha1 <- "85572c8d56504d8bba2afb32d7d7df35fb127ab8"

  expect_true(MALDIquantForeign:::.testChecksum(f, md5, algo="md5"))
  expect_true(MALDIquantForeign:::.testChecksum(f, sha1, algo="sha1"))
  expect_false(suppressWarnings(
                 MALDIquantForeign:::.testChecksum(f, "12345", algo="sha1")))
  expect_warning(MALDIquantForeign:::.testChecksum(f, "12345", algo="sha1"),
                 "Stored and calculated sha1 sums do not match")
  expect_warning(MALDIquantForeign:::.testChecksum(f, "12345", algo="md5"),
                 "Stored and calculated md5 sums do not match")
})
