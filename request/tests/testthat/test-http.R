context("http")

test_that("http works", {
  skip_on_cran()

  aa1 <- api("http://api.plos.org/search") %>%
    api_query(q = ecology, wt = json, fl = 'id,journal')

  aa2 <- api("http://api.plos.org/search") %>%
    api_query(q = ecology, wt = json, fl = 'id,journal') %>%
    peep()

  expect_is(aa1, "list")
  expect_is(aa2, "req")
  expect_identical(aa1, aa2 %>% http)

  expect_identical(http(api("https://api.github.com")),
    api("https://api.github.com") %>% http
  )

  x <- api("http://httpbin.org/post") %>%
      api_body(x = "A simple text string") %>%
      http("POST")
  expect_is(x, "list")
})

test_that("http fails well", {
  skip_on_cran()

  expect_error(http(), "argument \"req\" is missing")
  expect_error(http(api("https://api.github.com"), method = "FART"),
               "method must be one of GET or POST")
})
