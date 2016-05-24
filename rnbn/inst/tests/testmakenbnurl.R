context("Test makenbnurl")
test_that("makenbnurl recognises services", {
    expect_true(makenbnurl(service="obs", tvks="NBNSYS0000007073") ==
                  "https://data.nbn.org.uk/api/taxonObservations?ptvk=NBNSYS0000007073")
    expect_output(makenbnurl(service="feature", feature="284443"), 
                  "https://data.nbn.org.uk/api/features/284443")
    expect_output(makenbnurl(service="f", feature="284443"), 
                  "https://data.nbn.org.uk/api/features/284443")
    expect_output(makenbnurl(service="FEAT", feature="284443"), 
                  "https://data.nbn.org.uk/api/features/284443")
    expect_output(makenbnurl(service="taxon", tvks="NBNSYS0000007073"), 
                  "https://data.nbn.org.uk/api/taxa/NBNSYS0000007073")
    expect_error(makenbnurl(tvks="NBNSYS0000007073"), "no service specified")
    expect_error(makenbnurl(service="x", tvks="NBNSYS0000007073"), "service not recognised")
})

test_that("makenbnurl checks for minimum parameters", {
    expect_error(makenbnurl(service="obs"), "One of tvks or gridRef is required")
    expect_error(makenbnurl(service="fea"), "feature parameter is required")
    expect_error(makenbnurl(service="tax"), "tvks parameter is required")
})

test_that("makenbnurl doesn't mind type of numeric params", {
    expect_that(makenbnurl(service="obs", tvks="NBNSYS0000007073", startYear=1990), 
                prints_text("https://data.nbn.org.uk/api/taxonObservations\\?ptvk=NBNSYS0000007073&startYear=1990"))
    expect_that(makenbnurl(service="obs", tvks="NBNSYS0000007073", endYear=2010), 
                prints_text("https://data.nbn.org.uk/api/taxonObservations\\?ptvk=NBNSYS0000007073&endYear=2010"))
    expect_output(makenbnurl(service="f", feature=284443), 
                  "https://data.nbn.org.uk/api/features/284443")
})

test_that("makenbnurl checks for malformed parameters", {
    expect_error(makenbnurl(service="obs", tvks="NBNSYS00000070"),"tvks parameter is incorrect")
    expect_error(makenbnurl(service="tax", tvks="Pelecocera tricincta"),"tvks parameter is incorrect")
    expect_error(makenbnurl(service="obs", tvks="NBNSYS00 00007073"),"tvks parameter is incorrect")
    expect_error(makenbnurl(service="feature", feature="Holme Fen"), "feature parameter is incorrect")
    expect_error(makenbnurl(service="feature", feature="NBNSYS000"), "feature parameter is incorrect")
    expect_error(makenbnurl(service="obs", tvks="NBNSYS0000007073", datasets="SGB"), "datasets parameter is incorrect")
    expect_error(makenbnurl(service="obs", tvks="NBNSYS0000007073", startYear="26"), "startYear parameter is incorrect")
    expect_error(makenbnurl(service="obs", tvks="NBNSYS0000007073", endYear="192013"), "endYear parameter is incorrect")     
    expect_error(makenbnurl(service="obs", tvks="NBNSYS0000007073", startYear="2000", endYear="1990"), "startYear cannot be later than endYear")
    expect_error(makenbnurl(service="obs", tvks=c("NBNSYS0000007094", "NBNSYS00172195")),"tvks parameter is incorrect")
    expect_error(makenbnurl(service="obs", tvks="NBNSYS0000007073", datasets=c("SGB00001","GA000483","GA000152","GA00306")), "datasets parameter is incorrect")
})

test_that("makenbnurl copes with variations for occurrences service", {
    expect_that(makenbnurl(service="obs", tvks="NBNSYS0000007073", datasets="SGB00001"), 
                prints_text("https://data.nbn.org.uk/api/taxonObservations\\?ptvk=NBNSYS0000007073&datasetKey=SGB00001"))
    expect_that(makenbnurl(service="obs", tvks="NBNSYS0000007073", startYear="1990"), 
                prints_text("https://data.nbn.org.uk/api/taxonObservations\\?ptvk=NBNSYS0000007073&startYear=1990"))
    expect_that(makenbnurl(service="obs", tvks="NBNSYS0000007073", endYear="2010"), 
                prints_text("https://data.nbn.org.uk/api/taxonObservations\\?ptvk=NBNSYS0000007073&endYear=2010"))
    expect_that(makenbnurl(service="obs", tvks="NBNSYS0000007073", datasets="SGB00001", startYear="1990", endYear="2010"), 
                prints_text("https://data.nbn.org.uk/api/taxonObservations\\?ptvk=NBNSYS0000007073&datasetKey=SGB00001&startYear=1990&endYear=2010"))
    expect_that(makenbnurl(service="obs", tvks=c("NBNSYS0000007094", "NBNSYS0000172195"), datasets="SGB00001", startYear="1990", endYear="2010"), 
                prints_text("https://data.nbn.org.uk/api/taxonObservations\\?ptvk=NBNSYS0000007094&ptvk=NBNSYS0000172195&datasetKey=SGB00001&startYear=1990&endYear=2010"))
    expect_that(makenbnurl(service="obs", tvks="NBNSYS0000007073", datasets=c("SGB00001","GA000483","GA000152","GA000306")), 
                prints_text("https://data.nbn.org.uk/api/taxonObservations\\?ptvk=NBNSYS0000007073&datasetKey=SGB00001&datasetKey=GA000483&datasetKey=GA000152&datasetKey=GA000306"))
    expect_that(makenbnurl(service="obs", tvks=c("NBNSYS0000007094", "NBNSYS0000172195"), datasets=c("SGB00001","GA000483","GA000152","GA000306"), startYear="1990", endYear="2010"), 
                prints_text("https://data.nbn.org.uk/api/taxonObservations\\?ptvk=NBNSYS0000007094&ptvk=NBNSYS0000172195&datasetKey=SGB00001&datasetKey=GA000483&datasetKey=GA000152&datasetKey=GA000306&startYear=1990&endYear=2010"))
})
    