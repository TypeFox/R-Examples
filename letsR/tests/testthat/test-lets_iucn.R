context("Test for lets.iucn")

sp <- c("Musonycteris harrisoni", "Ailuropoda melanoleuca",
        "Cebus flavius")


test_that("lets.iucn works fine, one species", {
  skip_on_cran()
  testiucn <- lets.iucn("Panthera tigris")
  expect_equal(class(testiucn), "data.frame")
  testiucn2 <- lets.iucn("Panthera onca")
  expect_equal(class(testiucn2), "data.frame")
})

test_that("lets.iucn works fine, one species, count = TRUE", {
  skip_on_cran()
  testiucn <- lets.iucn("Panthera tigris", count = TRUE)
  expect_equal(class(testiucn), "data.frame")
  testiucn2 <- lets.iucn("Panthera onca", count = TRUE)
  expect_equal(class(testiucn2), "data.frame")
})


test_that("lets.iucn works fine, multiple species", {
  skip_on_cran()
  testiucn <- lets.iucn(sp)
  expect_equal(class(testiucn), "data.frame")
  testiucn2 <- lets.iucn(sp, count = TRUE)
  expect_equal(class(testiucn2), "data.frame")
  expect_equal(testiucn, testiucn2)
})
