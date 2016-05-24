context("exportTab")

m <- createMassSpectrum(mass=1:5, intensity=6:10)
r <- data.frame(V1=1:5, V2=6:10)

test_that("exportTab", {
  temp <- tempfile()
  MALDIquantForeign:::.exportTab(m, file=temp)
  ## didn't work on win-builder.r-project.org
  ## (but on local linux and windows install (both R 2.15.2)
  #expect_equivalent(tools::md5sum(temp), tools::md5sum(file.path("data",
  #                                                               "ascii.txt")))
  expect_equal(read.table(temp), r)
})

test_that("exportCsv", {
  temp <- tempfile()
  MALDIquantForeign:::.exportCsv(m, file=temp)
  ## didn't work on win-builder.r-project.org
  ## (but on local linux and windows install (both R 2.15.2)
  #expect_equivalent(tools::md5sum(temp), tools::md5sum(file.path("data",
  #                                                               "csv1.csv")))
  colnames(r) <- c("mass", "intensity")
  expect_equal(read.csv(temp), r)
})
