context ("Remote API Client: avisApiBusOrden")

response <- .avisApiBusOrden()

test_that(".avisApiBusOrden returns data", { 

	expect_true (length(response) > 0)
})

test_that(".avisApiBusOrden returns double", { 

	expect_true (is.double(response))
})
