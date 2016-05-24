context("Test for lets.iucn.ha")

sp <- c("Musonycteris harrisoni", "Ailuropoda melanoleuca",
        "Cebus flavius")


test_that("lets.iucn.ha works fine, one species", {
  skip_on_cran()
  testiucn <- lets.iucn.ha("Panthera tigris")
  expect_equal(class(testiucn), "data.frame")
  testiucn2 <- lets.iucn.ha("Panthera onca")
  expect_equal(class(testiucn2), "data.frame")
})

test_that("lets.iucn.ha works fine, one species, count = TRUE", {
  skip_on_cran()
  testiucn <- lets.iucn.ha("Panthera tigris", count = TRUE)
  expect_equal(class(testiucn), "data.frame")
  testiucn2 <- lets.iucn.ha("Panthera onca", count = TRUE)
  expect_equal(class(testiucn2), "data.frame")
})


test_that("lets.iucn.ha works fine, multiple species", {
  skip_on_cran()
  testiucn <- lets.iucn.ha(sp)
  expect_equal(class(testiucn), "data.frame")
  testiucn2 <- lets.iucn.ha(sp, count = TRUE)
  expect_equal(class(testiucn2), "data.frame")
  expect_equal(testiucn, testiucn2)
})
