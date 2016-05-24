
context('GetComps')

test_that("'getURL' errors are handled gracefully", {
    set_zillow_web_service_id('ZWSID')

    with_mock(
        getURL = function(...) {stop('Cryptic getURL error')},
        expect_error(GetComps(zpid = 48749425, count = 5), "Zillow API call with request '.+' failed with Error in RCurl::getURL\\(request\\): Cryptic getURL error"),
        .env = 'RCurl'
    )
})
