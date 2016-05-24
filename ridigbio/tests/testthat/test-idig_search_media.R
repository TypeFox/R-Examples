##' These tests are more media-specific, for fuller searching tests, see 
##' test-idig_search_records.R
##' 
context("test idig_search_media")

genus <- "cortinarius"
rq <-list("genus"=genus)
mq <- list("dqs"=list("type"="range", "gte"=0.2, "lte"=0.8))
fields <- c('uuid', 'dqs', 'hasSpecimen', 'data.ac:accessURI')
u <- "000113ce-84d4-467c-9fe3-13191596865e"

test_that("full results for rq searches return", {
  testthat::skip_on_cran()
  df <- idig_search_media(rq=rq, limit=6000)
  
  expect_that(df, is_a("data.frame"))
  expect_that(nrow(df) > 5000, is_true())
  expect_that(which(df$uuid == u) > 0, 
              is_true())
})

test_that("full results for rq & mq queries together return",{
  testthat::skip_on_cran()
  df <- idig_search_media(rq=rq, mq=mq, limit=6000)
  
  expect_that(df, is_a("data.frame"))
  expect_that(nrow(df) > 5000, is_true())
  expect_that(which(df$uuid == u) > 0, 
              is_true())
})

test_that("full results for mq queries return", {
  testthat::skip_on_cran()
  df <- idig_search_media(mq=mq, limit=6000)
  
  expect_that(df, is_a("data.frame"))
  expect_that(nrow(df) > 5000, is_true())
  expect_that(which(df$uuid == u) > 0, 
              is_true())
  expect_that(attributes(df)[["itemCount"]] > 5000, is_true())
  expect_that(length(attributes(df)[["attribution"]]) > 2 , is_true())
})

test_that("limits and custom fields return", {
  testthat::skip_on_cran()
  df <- idig_search_media(mq=mq, fields=fields, limit=10)
  
  expect_that(nrow(df) == 10, is_true())
  expect_that(ncol(df) == length(fields), is_true())
  expect_true(any(!is.null(df[["uuid"]])) && 
              any(df[["uuid"]] != "NA"))
  expect_true(any(!is.null(df[["data.ac:accessURI"]])) && 
              any(df[["data.ac:accessURI"]] != "NA"))
})
