# tests for avisHasSpecies in rAvis

context("avisHasSpecies")

test_that("Pica pica is in the database, Pica pic is not",{ 
	
    nameraw<- "Pica pica"
    expect_true(avisHasSpecies (nameraw))

    nameraw<- "Pica pic"
    expect_false(avisHasSpecies (nameraw))
})
