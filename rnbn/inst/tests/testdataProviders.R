context("Test dataProviders")

# login
load('~/rnbn_test.rdata')
nbnLogin(username = UN_PWD$username, password = UN_PWD$password)

test_that("Errors given", {
    expect_error(dataProviders(), 'datasets parameter cannot be NULL') 
})

test_that("Warnings given", {
    expect_warning(dataProviders(datasets='foo'), 'Dataset foo returned no taxa, check this is spelt correctly') 
})

test_that("Providers are returned", {
    expect_is(dps <- dataProviders(datasets=c('GA000312','GA000426')), 'data.frame')
    expect_that(nrow(dps) == 2, is_true())
    expect_that('Dr David Roy' %in% dps$contactName, is_true())
})
