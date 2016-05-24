context("avisSpeciesId")

# TODO: base on mocked .avisApiBusOrden

test_that("Correct species name is loaded", { 

  	nameraw <- "Pica pica"

  	expectedId <- 480

  	expect_is(avisSpeciesId (nameraw), "integer")
  	expect_identical (avisSpeciesId (nameraw), as.integer (expectedId))
})

test_that("Wrong species name throws an error", { 

	expect_error (avisSpeciesId ("Unexistent especies"))
})
