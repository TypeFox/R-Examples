context("Test getFeature")

# login
load('~/rnbn_test.rdata')
nbnLogin(username = UN_PWD$username, password = UN_PWD$password)

test_that("Get details for a feature", {
    expect_equal(as.character(getFeature("97479")['label']), "SN413499")    
})
