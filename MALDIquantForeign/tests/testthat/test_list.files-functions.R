context("list.files")


test_that("list.files-functions", {
  path <- normalizePath(system.file("exampledata", package="MALDIquantForeign"))
  expect_identical(MALDIquantForeign:::.list.files(path, pattern="tiny1-c"),
                   normalizePath(file.path(path,
                                           c("tiny1-centroided.mzML1.1.mzML",
                                             "tiny1-centroided.mzXML3.0.mzXML",
                                             "tiny1-compressed.mzML1.1.mzML",
                                             "tiny1-compressed.mzXML3.0.mzXML"))))
  expect_identical(MALDIquantForeign:::.list.files(path, pattern="tiny1-c",
                                                   excludePattern="\\.mzXML$"),
                   normalizePath(file.path(path,
                                           c("tiny1-centroided.mzML1.1.mzML",
                                             "tiny1-compressed.mzML1.1.mzML"))))
})
