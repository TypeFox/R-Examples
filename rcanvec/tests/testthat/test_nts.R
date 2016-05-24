#tests of NTS

context("Tests of NTS functions")

test_that("the nts() function retruns in the proper format", {
  ref <- nts('21h')
  ref2 <- nts('21h1')
  expect_equal(ref[1], "021", ref2[1])
  expect_equal(ref[2], "H", ref2[2])
  expect_equal(ref2[3], "01")
  expect_equal(length(ref), 2)
  expect_equal(length(ref2), 3)
})

test_that("multiple nts refs are output as a list", {
  refs <- nts('21h1', '21a16', '21A15')
  expect_is(refs, "list")
  expect_equal(length(refs), 3)
})

test_that("lat/lon locations locate the correct sheet", {
  expect_equal(nts(lat=45.2, lon=-64.32), nts("21h1"))
  expect_equal(nts(lat=c(45.2, 46.2), lon=c(-64.32, -64.81)), nts("21h1", "21i2"))
})

test_that("bounding boxes return correct sheets", {
  library(prettymapr)
  sheets <- nts(bbox=makebbox(45.125, -64.25, 44.875, -64.75))
  expect_equal(sheets, nts("021A15", "021H02", "021A16", "021H01"))
})

