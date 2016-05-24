context("vertmap")

test_that("vertmap works with vertsearch", {
  skip_on_cran()
  
  out <- vertsearch("mustela nigripes", verbose = FALSE, limit = 50)
  map1 <- vertmap(input = out, mapdatabase = "state")
  
  expect_is(map1, "gg")
  expect_is(map1$data, "data.frame")
  expect_equal(map1$labels$colour, "scientificname")
  expect_equal(as.character(map1$mapping$x), "long")
  expect_is(map1$mapping$x, "name")
})

test_that("vertmap works for maps not distinguished by color", {
  skip_on_cran()
  
  out <- vertsearch("mustela nigripes", verbose = FALSE, limit = 50)
  out$data$scientificname <- NULL
  map1 <- vertmap(input = out, mapdatabase = "state")
  
  expect_is(map1, "gg")
  expect_is(map1$data, "data.frame")
  expect_null(map1$labels$colour) # the difference here
  expect_equal(as.character(map1$mapping$x), "long")
  expect_is(map1$mapping$x, "name")
})

test_that("vertmap fails well", {
  skip_on_cran()
  
  expect_error(vertmap(), "Input must be of class list or data.frame")
  expect_error(vertmap(5), "Input must be of class list or data.frame")
})
