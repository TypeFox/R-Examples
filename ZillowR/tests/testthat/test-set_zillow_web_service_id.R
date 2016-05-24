
context('set_zillow_web_service_id')

test_that("'set_zillow_web_service_id' sets ZillowR-zws_id option", {
    zws_id <- getOption('ZillowR-zws_id')
    set_zillow_web_service_id('ZWSID')
    expect_equivalent(getOption('ZillowR-zws_id'), 'ZWSID')
    options('ZillowR-zws_id' = zws_id)
    rm(zws_id)
})
