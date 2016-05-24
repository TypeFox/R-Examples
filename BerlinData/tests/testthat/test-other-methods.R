context("other_methods")

test_that("download.CSV correctly loads and parses CSV", {
  csv.url <- './data/test-data.csv'
  csv.data <- download.CSV(csv.url)
  expect_error(download.CSV())
  expect_message(download.CSV(4))
  expect_warning(download.CSV('./data/test-data.json'))
  expect_equivalent(class(csv.data), "data.frame")
  expect_equivalent(dim(csv.data), c(4, 22))
})

test_that("download.JSON correctly loads and parses JSON", {
  json.url <- './data/test-data.json'
  json.data <- download.JSON(json.url)
  expect_error(download.JSON())
  expect_message(download.JSON(4))
  expect_message(download.JSON('./data/test-rss-feed.xml'))
  expect_equivalent(class(json.data), "data.frame")
  expect_equivalent(dim(json.data), c(4, 22))
})

test_that("download.XML correctly loads and parses XML", {
  xml.url <- './data/test-data.xml'
  xml.data <- download.XML(xml.url)
  expect_error(download.XML())
  expect_message(download.XML(4), 'Failed to download resource')
  expect_null(download.XML('./data/test-data.json'))
  expect_equivalent(class(xml.data), "data.frame")
  expect_equivalent(dim(xml.data), c(4, 22))
})

test_that("download.TXT correctly loads and parses TXT", {
  TXT.url <- './data/test-data.txt'
  TXT.data <- download.TXT(TXT.url)
  expect_error(download.TXT())
  expect_message(download.TXT(4), 'Failed to download resource')
  expect_warning(download.TXT('./data/test-data.json'))
  expect_equivalent(class(TXT.data), "data.frame")
  expect_equivalent(dim(TXT.data), c(4, 4))
})

test_that("download methods respond to verbosity options", {
  expect_that(download.default('foo', message.on.fail=FALSE), not(shows_message()))
  expect_that(download.CSV('foo', message.on.fail=FALSE), not(shows_message()))
  expect_that(download.JSON('foo', message.on.fail=FALSE), not(shows_message()))
  expect_that(download.XML('foo', message.on.fail=FALSE), not(shows_message()))
  expect_that(download.TXT('foo', message.on.fail=FALSE), not(shows_message()))
})
