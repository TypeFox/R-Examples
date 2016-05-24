context("UNFv6: Dates")
test_that("Partial dates (year-only) supported", {})
test_that("Partial dates (year-month) supported", {})

context("UNFv6: Datetimes")
test_that("Examples from v6 specification",{
    expect_equal(unf6("2014-08-22T16:51:05Z"), 
                 unf6(strptime("2014-08-22T16:51:05Z", "%Y-%m-%dT%H:%M:%OSZ", tz="UTC"), timezone="UTC"))
    #expect_equal(unf6("2012-06-10T14:29:00"), 
    #             unf6(strptime("2012-06-10T14:29:00", "%Y-%m-%dT%H:%M:%OS", tz="UTC"), timezone=""))
})

test_that("UNFs differ by timezone", {
    expect_false(identical(unf6(strptime("2014-08-22T16:51:05Z", "%Y-%m-%dT%H:%M:%OSZ", tz="UTC"), timezone="UTC"), 
                 unf6(strptime("2014-08-22T16:51:05Z", "%Y-%m-%dT%H:%M:%OSZ", tz="UTC"), timezone="US/Eastern")))
    expect_false(identical(unf6(strptime("2014-08-22T16:51:05Z", "%Y-%m-%dT%H:%M:%OSZ", tz="UTC"), timezone="UTC"), 
                 unf6(strptime("2014-08-22T16:51:05Z", "%Y-%m-%dT%H:%M:%OSZ", tz="US/Eastern"), timezone="UTC")))
    
})

test_that("Correct UNF for UTC timezones", {})

test_that("Tests of `decimal_seconds` rounding parameter", {
    expect_equal(unf6(as.POSIXct(1408726265.12345, origin="1970-01-01"), decimal_seconds=0),
                 unf6(as.POSIXct(1408726265, origin="1970-01-01")))
    expect_equal(unf6(as.POSIXct(1408726265.00000, origin="1970-01-01")),
                 unf6(as.POSIXct(1408726265, origin="1970-01-01")))
    expect_equal(unf6(as.POSIXct(1408726265.00000123456, origin="1970-01-01")),
                 unf6(as.POSIXct(1408726265, origin="1970-01-01")))
    expect_equal(unf6(as.POSIXct(1408726265.054321, origin="1970-01-01"), decimal_seconds=1),
                 unf6(as.POSIXct(1408726265.012345, origin="1970-01-01"), decimal_seconds=1))
    expect_false(identical(unf6(as.POSIXct(1408726265.12345, origin="1970-01-01"), decimal_seconds=0),
                           unf6(as.POSIXct(1408726265.12345, origin="1970-01-01"), decimal_seconds=1)))
})
