context ("avisContributorsSummary")

test_that("avisContributorsSummary returns expected header",{

	response <- avisContributorsSummary()

	expectedNames <- c("UserId", "Observations","Species","Provinces","UTMs","Periods")

  	expect_equal (colnames (response), expectedNames)
})
