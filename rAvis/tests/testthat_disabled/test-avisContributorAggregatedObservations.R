context ("avisContributorAggregatedObservations")

test_that("Returns expected header", {

	response <- avisContributorAggregatedObservations(100)

	expectedNames <- c("SpeciesId", "Observations", "Number", "UTM.10x10", "Birdwatchers")

	expect_equal (names (response), expectedNames)
})
