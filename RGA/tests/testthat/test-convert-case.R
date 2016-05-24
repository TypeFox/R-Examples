context("Convert camelCase character to separated")
 
test <- c("CamelCase", "CamelCamelCase", "Camel2Camel2Case", "getHTTPResponseCode", "get2HTTPResponseCode", "HTTPResponseCode", "HTTPResponseCodeXYZ")
result <- c("camel.case", "camel.camel.case", "camel2.camel2.case", "get.http.response.code", "get2.http.response.code", "http.response.code", "http.response.code.xyz")

test_that("Identical results", {
    expect_identical(to_separated(test), result)
})
