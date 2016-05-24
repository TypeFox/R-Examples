context("compression")

z <- c("bz2", "bzip2", "gz", "lzma", "xz")
e <- c("zip", z, paste("tar", z, sep="."), "tar", "txt", "mzML")
f <- paste(letters[1:13], e, sep=".")

test_that(".isCompressed", {
  expect_identical(MALDIquantForeign:::.isCompressed(f),
                   c(rep(TRUE, 11), rep(FALSE, 3)))
})

test_that(".isTar", {
  expect_identical(MALDIquantForeign:::.isTar(f),
                   c(rep(FALSE, 6), rep(TRUE, 6), rep(FALSE, 2)))
})

test_that(".isZip", {
  expect_identical(MALDIquantForeign:::.isZip(f), c(TRUE, rep(FALSE, 13)))
})

test_that(".isPackedOrCompressed", {
  expect_identical(MALDIquantForeign:::.isPackedOrCompressed(f),
                   c(rep(TRUE, 12), rep(FALSE, 2)))
})

test_that(".uncompress supports single file compression by gunzip", {
  u <- MALDIquantForeign:::.uncompress(
    system.file(file.path("exampledata", "compressed", "csv1.csv.gz"),
                package="MALDIquantForeign"))
  f <- system.file(file.path("exampledata", "csv1.csv"),
                   package="MALDIquantForeign")
  expect_identical(readLines(u), readLines(f))
  expect_identical(MALDIquantForeign:::.uncompress("foobar.txt"), "foobar.txt")
  expect_error(MALDIquantForeign:::.uncompress("foobar.gz"),
               ".*foobar.gz.* doesn't exist!")
})

test_that(".uncompress supports tar compression by untar", {
  u <- list.files(MALDIquantForeign:::.uncompress(
                    system.file(
                      file.path("exampledata", "compressed", "csv.tar.gz"),
                      package="MALDIquantForeign")),
                  recursive=TRUE, pattern="^.*\\.csv$", full.names=TRUE)[1]
  f <- system.file(file.path("exampledata", "csv1.csv"),
                   package="MALDIquantForeign")
  expect_identical(readLines(u), readLines(f))
})

test_that(".uncompress supports zip compression by unzip", {
  u <- list.files(MALDIquantForeign:::.uncompress(
                    system.file(
                      file.path("exampledata", "compressed", "csv.zip"),
                      package="MALDIquantForeign")),
                  recursive=TRUE, pattern="^.*\\.csv$", full.names=TRUE)[1]
  f <- system.file(file.path("exampledata", "csv1.csv"),
                   package="MALDIquantForeign")
  expect_identical(readLines(u), readLines(f))
  expect_error(MALDIquantForeign:::.uncompress("foobar.zip"),
               "unzip failed!")
})

test_that(".cleanupUncompressedTmpFiles works", {
  n <- list.files(file.path(tempdir(), "MALDIquantForeign_uncompress"),
                  recursive=TRUE)
  expect_true(length(n) > 0)
  MALDIquantForeign:::.cleanupUncompressedTmpFiles()
  expect_false(file.exists(file.path(tempdir(),
                                     "MALDIquantForeign_uncompress")))
})

test_that("typical auto import", {
  f <- normalizePath(system.file(
    file.path("exampledata", "compressed", "csv1.csv.gz"),
                     package="MALDIquantForeign"))
  s <- createMassSpectrum(mass=1:5, intensity=6:10)
  i <- import(f)[[1]]
  metaData(i) <- list()
  expect_identical(i, s)
  expect_false(file.exists(file.path(tempdir(),
                                     "MALDIquantForeign_uncompress")))
})
