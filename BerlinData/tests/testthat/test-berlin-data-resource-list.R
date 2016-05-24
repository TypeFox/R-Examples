context("berlin_data_resource_list")

test_bdr <- structure(list(
  title = 'bar',
  url = './data/test-data.xml',
  format = 'XML'
), class="berlin_data_resource")
test_bdrl <- structure(list(
    test_bdr,
    test_bdr
  ),
  class="berlin_data_resource_list")

test_that("berlin-data-resource-list methods give correct output", {
  expect_true(is.berlin_data_resource_list(test_bdrl))
  expect_output(summary(test_bdrl), '2 resources')
  expect_equivalent(dim(as.data.frame(test_bdrl)), c(2, 5))
  expect_equivalent(class(test_bdrl[1]), "berlin_data_resource_list")
})

test_that("berlin-data-resource-list download works correctly", {
  test_data <- download(test_bdrl)
  expect_message(download(test_bdrl), 'Downloading 2 resources')
  expect_equivalent(class(test_data), "list")
  expect_equivalent(class(test_data[[1]]), "data.frame")
  test_bdrl[[1]]$url <- ''
  expect_message(download(test_bdrl), 'Failed to download resource')
  expect_message(download(test_bdrl), 'Downloaded 1 of 2 resources')
})
