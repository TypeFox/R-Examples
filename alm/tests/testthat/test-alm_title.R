# tests for alm_title fxn in alm
context("alm_title")

test_that("almtitle returns the correct value", {
	expect_match(suppressWarnings(alm_title(doi='10.1371/journal.pbio.0000012')), "Genome-Wide RNAi of")
	expect_match(suppressWarnings(alm_title(doi='10.1371/journal.pbio.1001357')), "Niche-Associated Activation")
})

test_that("alm_title returns the correct class", {
	expect_is(suppressWarnings(alm_title(doi='10.1371/journal.pbio.0000012')), "character")
	expect_is(suppressWarnings(alm_title(doi='10.1371/journal.pbio.1001357')), "character")
})

test_that("alm_title returns a title of the correct length", {
	expect_that(nchar(suppressWarnings(alm_title(doi='10.1371/journal.pbio.1001357'))), equals(108))
	expect_that(nchar(suppressWarnings(alm_title(doi='10.1371/journal.pbio.0000012'))), equals(111))
})
