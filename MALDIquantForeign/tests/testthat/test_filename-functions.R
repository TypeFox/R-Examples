context("filename")

test_that(".cleanFilename", {
  expect_identical(MALDIquantForeign:::.cleanFilename(
                   "/home/a:/\"foo&bar\"/g.\\23!/ foo-bar?.txt"),
                   "_home_a_foo_bar_g_23_foo_bar_txt")
})

test_that("file extension is returned", {
  expect_identical(MALDIquantForeign:::.fileExtension("~/foo.txt"), "txt")
  expect_identical(MALDIquantForeign:::.fileExtension(
                     c("/etc/a.conf", "b.pdf")), c("conf", "pdf"))
})

test_that("file name is returned", {
  expect_identical(MALDIquantForeign:::.fileExtension("~/foo"), "foo")
})

test_that("path without extension is returned", {
  expect_identical(MALDIquantForeign:::.withoutFileExtension(
                     c("~/foo", "/home/user/xyz.tar.gz", "/tmp/bar.txt")),
                   c("~", "/home/user/xyz", "/tmp/bar"))
})

test_that("file extension is changed", {
  expect_identical(MALDIquantForeign:::.changeFileExtension(
                      c("/home/user/xyz.tar.gz", "/tmp/bar.txt"),
                      c("txt", "csv")),
                   c("/home/user/xyz.txt", "/tmp/bar.csv"))
})

test_that(".cutFilenames", {
  expect_identical(MALDIquantForeign:::.cutFilenames(
                     c("/home/user/foo.bar", "/home/user/xyz.tar.gz")),
                   c("foo.bar", "xyz.tar.gz"))

  expect_identical(MALDIquantForeign:::.cutFilenames(
                     c("/home/user/foo.bar", "/home/user/foo.bar")),
                   c("foo.bar", "foo.bar"))
})

test_that(".composeFilenames", {
  s <- createMassSpectrum(mass=1:5, intensity=1:5,
                          metaData=list(file="/foo/bar.txt"))
  expect_identical(MALDIquantForeign:::.composeFilename(s), "/foo/bar.csv")
  expect_identical(MALDIquantForeign:::.composeFilename(s, fileExtension="xml"),
                   "/foo/bar.xml")
  metaData(s) <- list(fullName="foo")
  expect_identical(MALDIquantForeign:::.composeFilename(s), "foo.csv")
  metaData(s) <- list(fullName=c("foo", "bar"))
  expect_identical(MALDIquantForeign:::.composeFilename(s), "foo_bar.csv")
})

test_that(".uniqueBaseFilenames", {
  expect_identical(MALDIquantForeign:::.uniqueBaseFilenames(
                     c("/home/user/foo.bar", "/home/user/foo.bar"),
                     fileExtension="txt"), c("foo_1.txt", "foo_2.txt"))
})

test_that(".make.unique", {
  expect_equal(MALDIquantForeign:::.make.unique(LETTERS[1:5]), LETTERS[1:5])
  expect_equal(MALDIquantForeign:::.make.unique(rep(LETTERS[1:5], each=2)),
               paste(rep(LETTERS[1:5], each=2), 1:2, sep="_"))
  expect_equal(MALDIquantForeign:::.make.unique(rep(LETTERS[1:2], each=10)),
               sprintf("%s_%02d", rep(LETTERS[1:2], each=10), 1:10))
})
