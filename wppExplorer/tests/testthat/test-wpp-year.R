context('Setting WPP year')

test_that('WPP year changes correctly', {
	set.wpp.year(2010)
	data <- wpp.indicator('leF') # contains observed data and predictions
	expect_true(setequal(unique(data$Year), seq(1955, 2100, by=5)))
	expect_equal(length(unique(data$charcode)), 195) # 195 countries
	
	set.wpp.year(2012)
	expect_false('leF' %in% ls(wppExplorer:::wpp.data.env)) # leF should be deleted from the environment
	data <- wpp.indicator('leF') # contains observed data, predictions and aggregations
	expect_true(setequal(unique(data$Year), seq(1955, 2100, by=5)))
	expect_equal(length(unique(data$charcode)), 235) 
	
	set.wpp.year(2015)
	expect_false('mig' %in% ls(wppExplorer:::wpp.data.env)) # mig should not be in the environment
	data <- wpp.indicator('mig') 
	expect_true(setequal(unique(data$Year), seq(1955, 2100, by=5)))
	expect_equal(length(unique(data$charcode)), 235) 
})