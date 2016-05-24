# tests for bison_tax fxn in taxize
context("bison_tax")

test_that("bison_tax returns the correct ...", {
  skip_on_cran()
  
  out1 <- bison_tax(query="*bear")
  out3 <- bison_tax(query="black bear", exact=TRUE)
  out4 <- bison_tax(query="helianthus", method="scientificName")
  
  # values
  expect_that(out1$numFound, equals(12))
  expect_that(out1$facets, equals(NULL))
  expect_that(out3$names$vernacularName, equals("black bear"))
  expect_that(out4$facets, equals(NULL))
  
  # class
  expect_that(out1, is_a("list"))
  expect_that(out3, is_a("list"))
  expect_that(out3$names, is_a("data.frame"))
  expect_that(out4, is_a("list"))
})
