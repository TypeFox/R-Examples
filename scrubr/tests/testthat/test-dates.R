context("dates functions")

df <- sample_data_1

test_that("standardizing dates works", {
  skip_on_cran()

  aa <- dframe(df) %>% date_standardize() %>% .$date %>% .[1]
  bb <- dframe(df) %>% date_standardize("%Y/%m/%d") %>% .$date %>% .[1]
  cc <- dframe(df) %>% date_standardize("%d%b%Y") %>% .$date %>% .[1]
  dd <- dframe(df) %>% date_standardize("%Y") %>% .$date %>% .[1]
  ee <- dframe(df) %>% date_standardize("%y") %>% .$date %>% .[1]

  expect_is(aa, "character")
  expect_is(bb, "character")
  expect_is(cc, "character")
  expect_is(dd, "character")
  expect_is(ee, "character")

  expect_match(aa, "^[0-9]{4}-[0-9]{2}-[0-9]{2}$")
  expect_match(bb, "^[0-9]{4}/[0-9]{2}/[0-9]{2}$")
  expect_match(cc, "^[0-9]{2}[A-Za-z]{3}[0-9]{4}$")
  expect_match(dd, "^[0-9]{4}$")
  expect_match(ee, "^[0-9]{2}$")
})

test_that("dropping date rows without data works", {
  skip_on_cran()

  aa <- dframe(df) %>% date_missing()

  expect_is(df, "data.frame")
  expect_is(aa, "data.frame")

  expect_equal(NROW(df), 1500)
  expect_equal(NROW(aa), 1498)
})

test_that("creating date rows data works", {
  skip_on_cran()

  df <- sample_data_2
  aa <- dframe(df) %>% date_create(year, month, day)

  expect_is(df, "data.frame")
  expect_is(aa, "data.frame")

  expect_equal(NCOL(df), 7)
  expect_equal(NCOL(aa), 8)
  expect_false(any(grepl("date", names(df))))
  expect_that(any(grepl("date", names(aa))), is_true())
})
