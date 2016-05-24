context("API")

dummy_res <- structure(list(
  content = raw(0),
  url = '',
  headers = list(`Content-Type` = "application/json;charset=utf-8"),
  status_code = 200
),
class = "response")

test_that("estat_getStatsList processes the API response as expected", {
  with_mock(
    `httr::GET` = function(...)
      purrr::update_list(dummy_res, content = readRDS("content_getStatsList.rds")),
    expect_identical(
      estat_getStatsList(
        appId = "XXXX",
        searchWord = "CD",
        limit = 3
      ),
      readRDS("result_getStatsList.rds")
    )
  )
})

test_that("estat_getMetaInfo processes the API response as expected", {
  with_mock(
    `httr::GET` = function(...)
      purrr::update_list(dummy_res, content = readRDS("content_getMetaInfo.rds")),
    expect_identical(
      estat_getMetaInfo(appId = "XXXX", statsDataId = "0003065345"),
      readRDS("result_getMetaInfo.rds")
    )
  )
})

test_that("estat_getStatsData processes the API response as expected", {
  with_mock(
    `httr::GET` = function(...)
      purrr::update_list(dummy_res, content = readRDS("content_getStatsData.rds")),
    expect_identical(
      estat_getStatsData(
        appId = "XXXX",
        statsDataId = "0003065345",
        cdCat01 = c("008", "009", "010"),
        limit = 3
      ),
      readRDS("result_getStatsData.rds")
    )
  )
})
