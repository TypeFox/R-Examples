context("pg_data")

test_that("works well", {
  skip_on_cran()

  pg_cache_clear(prompt = FALSE)

  aa <- pg_data(doi = '10.1594/PANGAEA.807580')

  expect_is(aa, "list")
  expect_is(aa[[1]], "pangaea")
  expect_is(unclass(aa[[1]]), "list")
  expect_named(unclass(aa[[1]]), c('doi', 'citation', 'meta', 'data'))
  expect_is(unclass(aa[[1]])$data, "tbl_df")

  expect_equal(pg_cache_list(), '10_1594_PANGAEA_807580.txt')
})

test_that("fails well", {
  skip_on_cran()

  expect_error(pg_data(), "\"doi\" is missing")
  expect_error(pg_data(5), "not of right form")
})
