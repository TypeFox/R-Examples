context('Indicators')

set.year <- function(wpp.year) {
	if(!get.wpp.year()==wpp.year)
		set.wpp.year(wpp.year)
}

test_that('fertility age profile has the right dimension', {
	prof <- wppExplorer:::get.age.profile.fert(2000, c('FR', 'AF'))
	expect_true(all(dim(prof) == c(2*7, 4))) # 2 countries by 7 age groups
	prof <- wppExplorer:::get.age.profile.fert(2020, 'GE')
	expect_true(all(dim(prof) == c(7, 4)))
	expect_true(all(prof$charcode == 'GE'))
})

test_that('age-specific fertility rate has the right value', {
	ind <- wpp.by.year(wpp.by.country(wpp.indicator('fertage', age="20-24"), 'FR'), 2015)$value
	data(tfr, package="wpp2015")
	data(percentASFR, package="wpp2015")
	tfrFR <- subset(tfr, country_code==250)["2010-2015"]
	asfrFR <- subset(percentASFR, country_code==250 & age=="20-24")["2010-2015"]
	expect_true(ind == tfrFR*asfrFR/100)
})


test_that('mean age of women in childbearing age is put into the right format', {
	set.year(2012)	
	data <- wpp.by.country(wpp.indicator('meanageinchbearage'), 'US')
	expect_true(all(dim(data) == c(31, 2))) # 1 country for all years
	expect_true(setequal(data$Year, seq(1950, 2100, by=5)))
})

test_that('mortality values for high ages come out correctly', {
	set.year(2015)
	mx <- wpp.by.year(wpp.by.country(wpp.indicator('mortagesex', sex="M", age="100+"), 'FR'), 2015)$value
	expect_true(mx > 0.51)
	mx <- wpp.by.year(wpp.by.country(wpp.indicator('mortagesex', sex="M", age="100+"), 'FI'), 2015)$value
	expect_true(mx > 0.5 & mx < 0.51)
	mx2 <- wpp.by.year(wpp.by.country(wpp.indicator('mortagesex', sex="M", age="100"), 'FI'), 2015)$value
	expect_true(all.equal(mx, mx2))
	mx <- wpp.by.year(wpp.by.country(wpp.indicator('mortagesex', sex="F", age="110"), 'FI'), 2015)$value
	expect_true(mx > 1.2)
})

test_that('migration rate has the right value', {
	ind <- wpp.by.year(wpp.by.country(wpp.indicator('migrate'), 'AE'), 2025)$value
	data(popproj, package="wpp2015")
	data(migration, package="wpp2015")
	popUA <- sum(subset(popproj, country_code==784)[c("2020", "2025")])/2
	migUA <- subset(migration, country_code == 784)["2020-2025"]
	expect_true(ind == migUA*200/popUA)
})