context("API calling")

if (run.integration.tests) {
    test_that("Request headers", {
        r <- try(crGET("http://httpbin.org/gzip"))
        expect_true(r$gzipped)
        expect_true(grepl("gzip", r$headers[["Accept-Encoding"]]))
        expect_true(grepl("rcrunch", r$headers[["User-Agent"]]))
    })

    test_that("crunchUserAgent", {
        expect_true(grepl("rcrunch", crunchUserAgent()))
        expect_true(grepl("libcurl", crunchUserAgent()))
        expect_error(crunchUserAgent("anotherpackage/3.1.4"),
            NA)
        expect_true(grepl("anotherpackage", crunchUserAgent("anotherpackage")))
    })

    with(test.authentication,
        test_that("API root can be fetched", {
            expect_true(is.shojiObject(getAPIroot()))
        }))

    test_that("API calls throw an error if user is not authenticated", {
        logout()
        expect_error(getAPIroot(),
            "You are not authenticated. Please `login\\(\\)` and try again.")
    })

    test_that("Deprecated endpoints tell user to upgrade", {
        expect_error(crGET("http://httpbin.org/status/410"),
            "The API resource at http://httpbin.org/status/410 has moved permanently. Please upgrade crunch to the latest version.")
    })
}
