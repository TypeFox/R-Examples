context("Test getOccurrences")

# login
load('~/rnbn_test.rdata')
nbnLogin(username = UN_PWD$username, password = UN_PWD$password)

test_that("Errors given", {
    expect_error(getOccurrences(tvks="badger", silent=T), 'Error in makenbnurl*') 
    expect_error(getOccurrences(tvks="NBNSYS0000002987", datasets="G3", silent=T), 'Error in makenbnurl*')
    expect_error(getOccurrences(tvks="NBNSYS0000002987", datasets="GA000373", startYear="1992", endYear="1991", silent=T), 'Error in makenbnurl*')
    expect_error(getOccurrences(tvks="NBNSYS0000002987", group="reptile", silent=T), 'Error in getOccurrences*')
    expect_error(getOccurrences(silent=T), 'Error in getOccurrences*')
})
dt <- getOccurrences(tvks="NBNSYS0000002987", datasets="GA000373", 
                     startYear="1990", endYear="1991", acceptTandC=TRUE, silent = TRUE)
test_that("Check single TVK search", {
    expect_that(nrow(dt) > 0, is_true()) 
    expect_that(sum(which(dt$absence == TRUE)), equals(0)) ## no absences
    expect_that(length(unique(dt$taxonVersionKey)), equals(1))
    expect_that(unique(dt$taxonVersionKey), equals('NBNSYS0000002987'))
    expect_that(length(unique(dt$datasetKey)), equals(1))
    expect_that(unique(dt$datasetKey), equals('GA000373'))
})
rm(dt)
dt <- getOccurrences(tvks=c("NBNSYS0000002987","NHMSYS0001688296","NHMSYS0000080210"),
                     startYear="1990", endYear="1991", acceptTandC=TRUE, silent = TRUE)
test_that("Check multi TVK search", {
    expect_that(nrow(dt) > 0, is_true()) 
    expect_that(sum(which(dt$absence == TRUE)), equals(0)) ## no absences
    expect_that(length(unique(dt$pTaxonVersionKey)), equals(3))
    expect_that(sort(unique(dt$pTaxonVersionKey)), equals(c('NBNSYS0000002987','NHMSYS0000080210','NHMSYS0001688296')))
})
rm(dt)
dt <- getOccurrences(group="quillwort", startYear="1990", endYear="1992",
                     VC="Shetland (Zetland)", acceptTandC=TRUE, silent = TRUE)
test_that("Check group search", {
    expect_that(nrow(dt) > 0, is_true()) 
    expect_that(sum(which(dt$absence == TRUE)), equals(0)) ## no absences
    expect_that(length(unique(dt$pTaxonVersionKey)), equals(2))
    expect_that(sort(unique(dt$pTaxonVersionKey)), equals(c('NBNSYS0000002008','NBNSYS0000002009')))
})
rm(dt)
dt <- getOccurrences(gridRef='SP00',startYear=2008,endYear=2008,acceptTandC=TRUE, silent = TRUE)
hecs <- unique(do.call(rbind, lapply(dt$location, FUN=gridRef, format='sq10km'))[,3])
test_that("Check gridRef only search", {
    expect_that(nrow(dt) > 0, is_true()) 
    expect_that(sum(which(dt$absence == TRUE)), equals(0)) ## no absences
    expect_that(length(hecs), equals(1))
    expect_that(hecs[[1]], equals('SP00'))
})
test_that("Check providers are returned", {
    expect_true('providers' %in% names(attributes(dt)))
    expect_is(attr(dt,'providers'), 'data.frame')
    expect_true(ncol(attr(dt,'providers')) >= 7)
    expect_true('Biological Records Centre' %in% attr(dt,'providers')$name)
})
rm(dt)
WCresults <- getOccurrences(tvks = 'NHMSYS0000332741', startYear = 1999,
                            endYear = 1999, attributes = TRUE, acceptTandC = TRUE, silent = TRUE)
test_that("Attributes are returned as expected", {    
    expect_true('attributes.Comment' %in% names(WCresults))
    expect_true(length(names(WCresults)[grepl('attributes.',names(WCresults))]) >=3)
    expect_is(WCresults, 'data.frame')
})