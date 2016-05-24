context("berlin_data_dataset")

test_bdr <- structure(list(
  title = 'bar',
  url = './data/test-data.xml',
  format = 'XML'
), class="berlin_data_resource")
test_bdd <- structure(list(
  title = 'foo',
  resources = structure(list(
    test_bdr
  ), class="berlin_data_resource_list")
), class="berlin_data_dataset")

test_that("berlin_data_dataset methods give correct output", {
  expect_true(is.berlin_data_dataset(test_bdd))
  expect_output(summary(test_bdd), '1 resource')
  expect_equivalent(dim(as.data.frame(test_bdd)), c(1, 6))
  expect_equivalent(resources(test_bdd)[[1]], test_bdr)
  expect_equivalent(test_bdd, getDatasetMetaData(test_bdd))
})

test_that("berlin_data_dataset download method functions correctly", {
  test_data <- download(test_bdd)
  expect_message(download(test_bdd), 'Downloading all resources for dataset')
  expect_equivalent(class(test_data), "list")
  expect_equivalent(class(test_data[[1]]), "data.frame")
  expect_equivalent(dim(test_data[[1]]), c(4, 22))
  test_bdd$resources <- list()
  expect_message(download(test_bdd), 'No resources located')
  expect_that(download(test_bdd, message.on.succeed=FALSE), 
              not(shows_message('Downloaded XML')))
})
