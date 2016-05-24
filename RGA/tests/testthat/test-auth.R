context("Authorization")

check_token <- function() {
    if (!file.exists(".ga-token.rds"))
        skip("Access token not available")
    suppressMessages(authorize(cache = ".ga-token.rds"))
}

test_that("Validate access token", {
    check_token()
    expect_true(validate_token(get_token()))
})

test_that("Test API request with access token", {
    check_token()
    url <- "https://www.googleapis.com/analytics/v3/management/accounts?max-results=1"
    req <- httr::GET(url, httr::config(token = get_token()))
    expect_equal(httr::status_code(req), 200L)
})
