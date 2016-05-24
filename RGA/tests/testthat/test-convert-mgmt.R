context("Convert Management API response")

mgmt_data <- structure(list(
    kind = "analytics#profiles",
    totalResults = 8L,
    items = structure(list(
        id = c(83638429L, 83639328L, 84435797L, 90904152L, 91911821L),
        accountId = c(26921269L, 26921269L, 26921269L, 26921269L, 26921269L),
        webPropertyId = c("UA-26921269-2", "UA-26921269-3", "UA-26921269-4", "UA-26921269-5", "UA-26921269-6"),
        name = c("All data", "All data", "All data", "All data", "All data"),
        currency = c("USD", "RUB", "USD", "USD", "USD"),
        timezone = c("Europe/Moscow", "Europe/Moscow", "Europe/Moscow", "Europe/Moscow", "Europe/Moscow"),
        websiteUrl = c("http://example.com", "http://example.com", "http://example.com", "http://example.com", "http://example.com"),
        type = c("WEB", "WEB", "WEB", "WEB", "WEB"),
        created = c("2014-03-19T03:24:37.347Z", "2014-03-19T03:33:41.892Z", "2014-04-06T05:06:16.123Z", "2014-09-09T02:19:02.907Z", "2014-10-01T16:26:27.409Z"),
        updated = c("2014-12-25T16:18:25.789Z", "2014-12-25T16:17:52.015Z", "2014-12-25T16:18:25.789Z", "2014-12-25T16:18:25.789Z", "2014-12-25T16:18:25.789Z"),
        eCommerceTracking = c(FALSE, FALSE, FALSE, FALSE, FALSE),
        permissions.effective = list(
            c("COLLABORATE", "EDIT", "MANAGE_USERS", "READ_AND_ANALYZE"),
            c("COLLABORATE", "EDIT", "MANAGE_USERS", "READ_AND_ANALYZE"),
            c("COLLABORATE", "EDIT", "MANAGE_USERS", "READ_AND_ANALYZE"),
            c("COLLABORATE", "EDIT", "MANAGE_USERS", "READ_AND_ANALYZE"),
            c("COLLABORATE", "EDIT", "MANAGE_USERS", "READ_AND_ANALYZE"))),
            .Names = c("id", "accountId", "webPropertyId", "name", "currency", "timezone", "websiteUrl", "type", "created", "updated", "eCommerceTracking", "permissions.effective"),
            class = "data.frame",
            row.names = c(NA, 5L))),
    .Names = c("kind", "totalResults", "items"))

mgmt_df <- build_df(mgmt_data)

test_that("Result class", {
    expect_is(mgmt_df, "data.frame")
})

test_that("Data frame dimensions", {
    expect_equal(ncol(mgmt_df), 12L)
    expect_equal(nrow(mgmt_df), 5L)
})

test_that("Columns types", {
    expect_is(mgmt_df[, 1L], "integer")
    expect_is(mgmt_df[, 2L], "integer")
    expect_is(mgmt_df[, 11L], "logical")
})
