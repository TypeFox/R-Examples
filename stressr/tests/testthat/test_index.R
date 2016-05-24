# unit tests for stressr stress index data
require(testthat)
require(stressr)
require(lattice)

context("Stress Index Data")

# this setup mimics the network call function
# reading our abbreviated test file
x <- readHTMLTable("cfsi-data.xml",
                   header=TRUE,
                   skip.rows=1,
                   colClasses=c("character",rep("numeric",1)),
                   stringsAsFactors=FALSE)

cfsi <- x[[1]]
colnames(cfsi) <- c(
  "DATE",
  "CLEVELAND")

cfsi$DATE <- as.Date(cfsi$DATE,format="%m/%d/%Y")
cfsi <- xts(cfsi[,-1],order.by=cfsi[,1])
colnames(cfsi) <- "CLEVELAND"

idx <- list(df=cfsi,
           colors="#376e90",
           main="Financial Stress Index",
           ylab="(Z-Score)")
class(idx) <- "cfsi"


test_that("correctly loaded data frame", {
  expect_equal(nrow(idx$df),908)
  expect_equal(ncol(idx$df),1)
  expect_equal(colnames(idx$df),c("CLEVELAND"))
  expect_equal(as.numeric(last(idx$df)),-0.895)
  expect_true("xts" %in% class(idx$df))
})

test_that("correctly loaded list", {
  expect_equal(names(idx),c("df","colors","main","ylab"))
  expect_equal(length(idx$colors),1)
  expect_equal(idx$main,"Financial Stress Index")
  expect_equal(idx$ylab,"(Z-Score)")
})

test_that("performs xyplot",{
  expect_equal(class(xyplot(idx)),"trellis")
})

test_that("performs stress index plot",{
  expect_equal(class(stressIndexChart(idx)),"trellis")
  expect_equal(class(stressIndexChart(idx,range="2013")),"trellis")
  expect_equal(class(stressIndexChart(idx,showGradeRegions=FALSE)),"trellis")
  expect_equal(class(stressIndexChart(idx,range="2012",showGradeRegions=TRUE)),"trellis")
})
  