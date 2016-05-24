context("api_query")

test_that("api_query works", {
  skip_on_cran()

  expect_is(api_query(api("http://api.plos.org/search")), "req")

  aa <- api("http://api.plos.org/search") %>%
    api_query(q = ecology, wt = json, fl = 'id,journal') %>%
    peep

  bb <- api("http://api.plos.org/search") %>%
    api_query(q = ecology, wt = json, fl = id, fl = journal) %>%
    peep

  cc <- api("http://api.plos.org/search") %>%
    api_query_(q = "ecology", wt = "json", fl = 'id', fl = 'journal') %>%
    peep

  expect_is(aa, "req")
  expect_is(bb, "req")
  expect_is(cc, "req")

  expect_is(aa$url, "url")
  expect_is(bb$query, "list")

  expect_is(aa %>% http, "list")
  expect_is(bb %>% http, "list")
  expect_is(cc %>% http, "list")

  expect_identical(cc %>% http,
    api("http://api.plos.org/search") %>%
      api_query_(q = "ecology", wt = "json", fl = 'id', fl = 'journal')
  )
})

test_that("api_query fails well", {
  skip_on_cran()

  xx <- api("http://api.plos.org/search")

  expect_error(api_query(), "argument \".data\" is missing")
})
