context("berlin_data_resource")

test_that("download method works correctly", {
  test_bdr <- structure(list(
    url = './data/test-data.xml',
    format = 'XML'
    ), class="berlin_data_resource")
  data = download(test_bdr)
  expect_equivalent(class(data), "data.frame")
  expect_equivalent(dim(data), c(4, 22))
})

test_that("download method gives appropriate errors", {
  test_bdr <- structure(list(
    url = 'https://not.supported.url',
    format = 'UNSUPPORTED-FORMAT'
  ), class="berlin_data_resource")
  expect_message(download(test_bdr), 'does not support the https:// URL scheme')
  test_bdr$url <- './data/test-data.xml'
  expect_message(download.berlin_data_resource(test_bdr), 'does not currently support')
  expect_error(download.berlin_data_resource(4))
  expect_error(download.berlin_data_resource())
})

test_that("berlin_data_resource gives right output for methods", {
  test_bdr <- structure(list(
    url = './data/test-data.xml',
    format = 'XML'
  ), class="berlin_data_resource")
  expect_true(is.berlin_data_resource(test_bdr))
  expect_output(summary(test_bdr), 'Format: XML')
  expect_equivalent(dim(as.data.frame(test_bdr)), c(1, 3))
  expect_message(download(test_bdr, message.on.succeed=TRUE))
  expect_that(download(test_bdr, message.on.succeed=FALSE), not(shows_message()))
})
