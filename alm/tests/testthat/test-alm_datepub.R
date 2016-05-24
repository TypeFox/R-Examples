# tests for alm_datepub fxn in alm
context("alm_datepub")

test_that("alm_datepub returns the correct value", {
	expect_that(alm_datepub('10.1371/journal.pone.0026871'), matches("2012-03-26"))
	expect_that(alm_datepub('10.1371/journal.pone.0026871', get = 'year'), equals(2012))
})

test_that("alm_datepub returns the correct class", {
	expect_that(alm_datepub('10.1371/journal.pone.0026871'), is_a("character"))
	expect_that(alm_datepub('10.1371/journal.pone.0026871', get='year'), is_a("numeric"))
})
