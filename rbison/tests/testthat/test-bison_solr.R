# tests for bison_solr fxn in taxize
context("bison_solr")

test_that("bison returns the correct value", {
  skip_on_cran()
  
  out_1 <- bison_solr(scientificName='Ursus americanus', verbose = FALSE)
  out_2 <- bison_solr(scientificName='Ursus americanus', state_code='New Mexico', fl="scientificName", verbose = FALSE)
  out_3 <- bison_solr(scientificName='Ursus americanus', state_code='New Mexico', rows=50, fl="occurrence_date,scientificName", verbose = FALSE)
  out_5 <- bison_solr(scientificName='Helianthus annuus', rows=800, verbose = FALSE)
  out_5_map <- bisonmap(out_5)

  # values
  expect_that(out_1$facets$facet_queries, equals(NULL))
  expect_that(out_2$highlight, equals(NULL))
  expect_that(out_2$points[[1]][1], equals("Ursus americanus"))

  # class
  expect_is(out_1$points, "data.frame")
  expect_is(out_1, "bison_solr")
  expect_is(out_2, "bison_solr")
  expect_is(out_3, "bison_solr")
  expect_is(out_5_map, "list")
  expect_is(out_5_map$plot, "gg")
})
