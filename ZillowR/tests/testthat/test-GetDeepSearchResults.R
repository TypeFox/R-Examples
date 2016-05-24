
context('GetDeepSearchResults')

test_that("'getURL' errors are handled gracefully", {
    set_zillow_web_service_id('ZWSID')

    with_mock(
        getURL = function(...) {stop('Cryptic getURL error')},
        expect_error(GetDeepSearchResults(address = '2114 Bigelow Ave', citystatezip = 'Seattle, WA'), "Zillow API call with request '.+' failed with Error in RCurl::getURL\\(request\\): Cryptic getURL error"),
        .env = 'RCurl'
    )
})
