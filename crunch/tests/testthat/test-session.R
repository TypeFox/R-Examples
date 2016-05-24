context("Sessions")

with_mock_HTTP({
    test_that("session() returns a session object", {
        expect_true(is.list(session()))
        expect_true("datasets" %in% names(session()))
    })
})
