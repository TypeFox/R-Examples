context("berlin_data_dataset_info")

test_that("berlin_data_dataset_info gives correct output for methods", {
  dataset_url <- './data/data-datasetpage2.html'
  test_bddi <- structure(list(description = "foo", 
                              title = "bar", 
                              link = dataset_url), 
                         class = "berlin_data_dataset_info")
  expect_true(is.berlin_data_dataset_info(test_bddi))
  expect_message(download(test_bddi))
  expect_null(download(test_bddi))
  expect_equivalent(dim(as.data.frame(test_bddi)), c(1,3))
  expect_equivalent(getDatasetMetaData(test_bddi), getDatasetMetaData(test_bddi$link))
})