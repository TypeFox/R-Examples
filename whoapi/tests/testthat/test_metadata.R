context("Test retrieval of metadata")

test_that("Domain search ranking can be retrieved", {
  token <- whoapi_token("demokey")
  result <- try({domain_rank(token, "whoapi.com")}, silent = TRUE)
  if(!"try-error" %in% class(result)){
    expect_that(is.list(result), equals(TRUE))
    expect_that(names(result), equals(c("status","pr","alexa_reach","alexa_popularity","alexa_linksin", "email", "title")))
  }
})

test_that("domain search results can be retrieved", {
  token <- whoapi_token("demokey")
  result <- try({domain_search(token, "whoapi.com")}, silent = TRUE)
  if(!"try-error" %in% class(result)){
    expect_that(is.list(result), equals(TRUE))
    expect_that(names(result), equals(c("status","google_results","bing_results")))
  }
})

test_that("Metadata (in the HTML sense) can be retrieved", {
  token <- whoapi_token("demokey")
  result <- try({domain_metadata(token, "whoapi.com")}, silent = TRUE)
  if(!"try-error" %in% class(result)){
    expect_that(is.list(result), equals(TRUE))
    expect_that(names(result), equals(c("status","title","description")))
  }
})

test_that("Geolocation info can be retrieved", {
  token <- whoapi_token("demokey")
  result <- try({domain_location(token, "whoapi.com")}, silent = TRUE)
  if(!"try-error" %in% class(result)){
    expect_that(is.list(result), equals(TRUE))
    expect_that(names(result), equals(c("status","ip","geo_cc","geo_country","geo_city", "geo_latitude", "geo_longitude")))
  }
})