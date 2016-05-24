context("basic functionality")
test_that("the API works", {

  test_urls <- c("http://t.co/D4C7aWYIiA",
    "1.usa.gov/1J6GNoW",
    "ift.tt/1L2Llfr",
    "bit.ly/1GPr5w5",
    "http://l.dds.ec/1da152x",
    "http://l.rud.is/seven")

  expect_that(expand_urls(test_urls), is_a("data.frame"))
  expect_that(known_services(), is_a("data.frame"))

})
