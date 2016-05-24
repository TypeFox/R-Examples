context("flatten_query")

test_that("flatten_query works", {
  expect_equal(flatten_query(list(statsDataId = "0003103532", cdCat01 = c("010800130","010800140"))),
               list(statsDataId = "0003103532", cdCat01 = "010800130,010800140"))
})
