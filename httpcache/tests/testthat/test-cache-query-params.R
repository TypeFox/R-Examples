context("Cache with GET query parameters")

public({
    clearCache()
    test_that("Checking cache even with cache off doesn't fail on long query", {
        uncached({
            with_mock_HTTP({
                z <- GET("https://beta.crunch.io/api/users/",
                    query=list(query=rep("Q", 10000)))
            })
        })
        expect_true(is.numeric(z$response))
    })

    clearCache()
    test_that("cache gets set on GET even with long query", {
        with_mock_HTTP({
            GET("https://beta.crunch.io/api/users/",
                query=list(query=rep("Q", 10000)))
        })
        expect_identical(cacheKeys(),
            "https://beta.crunch.io/api/users/?HASHED_QUERY=38f0ed36c36e7c08ad375cc9a48d1364")
    })
    without_internet({
        test_that("Can read cache with query params even with no connection", {
            expect_identical(GET("https://beta.crunch.io/api/users/",
                query=list(query=rep("Q", 10000)))$response,
                33L)
        })
        test_that("Caching respects GET query parameters", {
            ## This is a cache miss because the query param is different
            expect_error(GET("https://beta.crunch.io/api/users/",
                query=list(a=1)))
        })
    })
})
