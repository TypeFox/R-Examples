context ("Remote API Client: avisApiBusEspecie")

response<- .avisApiBusEspecie()

test_that(".avisApiBusEspecie output is a dataframe",{ 
  
	expect_is(response, 'data.frame')
})

test_that(".avisApiBusEspecie returns expected header",{ 

	expectedHeader <- c("Observations", "Individuals", "UTM.10x10", "Birdwatchers")

	expect_equal (names(response), expectedHeader)
})

test_that(".avisApiBusEspecie returns records",{ 

	expect_equal (dim(o)[1] > 0)
})
