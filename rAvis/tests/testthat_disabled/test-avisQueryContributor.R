context ("avisQueryContributor")

# TODO: test based on mocked .avisApiFichaUsuario and .avisApiUsuarios

test_that("avisQueryContributor returns correct format",{ 

  	# get some valid contributor id
  	cs <- avisContributorsSummary()
  	id <- cs[[1]]

  	expect_is(avisQueryContributor(id), 'data.frame')        
})
