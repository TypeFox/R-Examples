# tests for bisonmap fxn in taxize
context("bisonmap")

test_that("bisonmap returns the correct class", {
  skip_on_cran()
  
  out <- bison(species="Aquila chrysaetos", count=10)
  map1 <- bisonmap(out)
  map2 <- bisonmap(out, tomap="county")
  map3 <- bisonmap(out, tomap="state")
  
  expect_is(map1, "gg")
  expect_is(map1$data, "data.frame")
  expect_is(map1$scales, "refClass")
  expect_is(map2, "gg")
  expect_is(map3, "gg")
})
