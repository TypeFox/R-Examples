context("Test datasetTaxa")

# login
load('~/rnbn_test.rdata')
nbnLogin(username = UN_PWD$username, password = UN_PWD$password)

test_that("Errors given", {
    expect_error(datasetTaxa(), 'datasets parameter cannot be NULL') 
})

test_that("Warnings given", {
    expect_warning(datasetTaxa(datasets='foo'), 'Dataset foo returned no taxa, check this is spelt correctly') 
})

test_that("Taxa are returned", {
    expect_is(taxa <- datasetTaxa(datasets=c('GA000312')), 'data.frame')
    expect_that(nrow(taxa) > 0, is_true())
    expect_that('Harmonia axyridis' %in% taxa$name, is_true())
})
