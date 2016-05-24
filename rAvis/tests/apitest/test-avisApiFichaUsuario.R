context ("Remote API Client: avisApiFichaUsuario")

testContributorId <- 370 # id of a contributor know to exist in remote db and with data

response <- .avisApiFichaUsuario (testContributorId) 

test_that(".avisApiFichaUsuario output is a dataframe",{ 

	expect_true(is.data.frame (response))
})

test_that(".avisApiFichaUsuario output has expected names",{ 

	# preliminar specification
	expectedNames <- c("SpeciesId", "Observations", "Number", "UTM.10x10", "Birdwatchers")

    expect_equal (names (response), expectedNames)
})

test_that(".avisApiFichaUsuario returns data",{ 

	# preliminar specification
	expectedNames <- c("SpeciesId", "Observations", "Number", "UTM.10x10", "Birdwatchers")

    expect_equal (names (response), expectedNames)
})
