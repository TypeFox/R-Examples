context("pg_search")

test_that("pg_search works and stuff, and stuff and things, also, it works", {
  skip_on_cran()

  aa <- pg_search(query='water')
  bb <- pg_search(query='water', count=2)
  cc <- pg_search(query='water', env="water")
  dd <- pg_search(query='water', mindate="2013-06-01", maxdate="2013-07-01")
  ee <- pg_search(query='water', bbox=c(-124.2, 41.8, -116.8, 46.1))
  ff <- pg_search(query='citation:Archer')

  expect_is(aa, "tbl_df")
  expect_named(aa, c('doi','score','size_datasets','citation','supplement_to'))
  expect_match(aa$doi, "10.1594")
  expect_is(aa$score, "numeric")
  expect_is(aa$size_datasets, "numeric")
  expect_is(aa$supplement_to, "character")

  expect_is(bb, "tbl_df")
  expect_is(cc, "tbl_df")
  expect_is(dd, "tbl_df")
  expect_is(ee, "tbl_df")
  expect_is(ff, "tbl_df")
})

test_that("fails well", {
  skip_on_cran()

  expect_error(pg_search(), "argument \"query\" is missing")
  expect_error(pg_search("water", count = "asdafd"), "count must be of class")
  expect_error(pg_search("water", env = 5), "env must be of class")
  expect_error(pg_search("water", mindate = 5), "mindate must be of class")
  expect_error(pg_search("water", maxdate = 5), "maxdate must be of class")
})
