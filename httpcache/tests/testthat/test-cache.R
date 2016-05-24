context("Caching")

test_that("Cache is on by default", {
    expect_true(is.null(getOption("httpcache.on")))
    expect_true(caching())
})
test_that("httpcache.on option affects caching", {
    on.exit(options(httpcache.on=NULL))
    options(httpcache.on=FALSE)
    expect_false(caching())
    options(httpcache.on=TRUE)
    expect_true(caching())
})

public({
    clearCache()

    test_that("Cache gets set on GET", {
        expect_identical(length(cacheKeys()), 0L)
        with_mock_HTTP({
            a <- GET("https://beta.crunch.io/api/datasets")
            b <- GET("https://beta.crunch.io/api/", query=list(user="me"))
        })
        expect_identical(length(cacheKeys()), 2L)
        expect_true("https://beta.crunch.io/api/datasets" %in% cacheKeys())
        expect_identical(a$response, 35L)
    })

    without_internet({
        test_that("When the cache is set, can read from it even with no connection", {
            ## Now read from cache
            expect_identical(GET("https://beta.crunch.io/api/datasets")$response,
                35L)
        })
        test_that("But uncached() prevents reading from the cache", {
            expect_error(uncached(GET("https://beta.crunch.io/api/datasets")),
                "GET https://beta.crunch.io/api/datasets")
        })
    })

    test_that("PUT busts cache", {
        ## Now bust cache
        with_mock_HTTP({
            expect_message(PUT("https://beta.crunch.io/api/datasets"),
                "PUT https://beta.crunch.io/api/datasets ")
        })
        ## See that it's no longer in the cache
        expect_identical(length(cacheKeys()), 1L)
        without_internet({
            expect_error(GET("https://beta.crunch.io/api/datasets"),
                "GET https://beta.crunch.io/api/datasets")
        })
    })

    test_that("PATCH busts cache", {
        without_internet({
            ## It's in the cache
            expect_identical(GET("https://beta.crunch.io/api/",
                query=list(user="me"))$response, 27L)
        })
        ## Now bust cache
        with_mock_HTTP({
            expect_message(PATCH("https://beta.crunch.io/api/"),
                "PATCH https://beta.crunch.io/api/ ")
        })
        ## See that it's no longer in the cache
        expect_identical(length(cacheKeys()), 0L)
        without_internet({
            expect_error(GET("https://beta.crunch.io/api/", query=list(user="me")),
                "GET https://beta.crunch.io/api/")
        })
    })

    test_that("cacheOff stops caching and clears existing cache", {
        with_mock_HTTP({
            GET("https://beta.crunch.io/api/datasets")
        })
        expect_identical(length(cacheKeys()), 1L)
        cacheOff()
        on.exit(cacheOn()) ## Turn it back on
        expect_identical(length(cacheKeys()), 0L)
        with_mock_HTTP({
            a <- GET("https://beta.crunch.io/api/datasets")
        })
        expect_identical(length(cacheKeys()), 0L)
        expect_identical(a$response, 35L)
    })
})
