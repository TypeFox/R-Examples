context("Test retrieval of direct information")

test_that("Test we can identify if a domain is taken", {
  token <- whoapi_token("demokey")
  result <- try({is_taken(token, "whoapi.com")}, silent = TRUE)
  if(!"try-error" %in% class(result)){
    expect_that(result, equals(TRUE))
  }

})

test_that("Test we can identify if a domain is blacklisted", {
  token <- whoapi_token("demokey")
  result <- try({is_blacklisted(token, "whoapi.com")}, silent = TRUE)
  if(!"try-error" %in% class(result)){
    expect_that(is.list(result), equals(TRUE))
    expect_that("blacklisted" %in% names(result), equals(TRUE))
  }
})

test_that("Test we can run whois against a domain", {
  token <- whoapi_token("demokey")
  result <- try({whois_info(token, "whoapi.com")}, silent = TRUE)
  if(!"try-error" %in% class(result)){
    expect_that(is.list(result), equals(TRUE))
    expect_that(c("contacts","whois_server","nameservers") %in% names(result), equals(c(TRUE,TRUE,TRUE)))
  }
})

test_that("Certificate information can be retrieved", {
  result <- try({certificate_info(token, "whoapi.com")}, silent = TRUE)
  if(!"try-error" %in% class(result)){
    expect_that(is.list(result), equals(TRUE))
    expect_that(c("status","certondomain") %in% names(result), equals(c(TRUE,TRUE)))
  }
})