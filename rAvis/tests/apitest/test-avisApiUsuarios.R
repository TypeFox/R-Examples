context ("Remote API Client: avisApiUsuarios")

response <- .avisApiUsuarios()

test_that(".avisApiUsuarios output is a matrix",{

  	expect_true(is.matrix (response))
})

test_that(".avisApiUsuarios returns expected header",{

	expectedNames <- c("UserId", "User","Observations","Species","Provinces","UTMs","Periods")

  	expect_equal (colnames (response), expectedNames)
})
