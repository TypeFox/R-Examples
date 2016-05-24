context("berlin-data-list")

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

test_bl1 <- structure(list(
  test_bdd,
  test_bdd,
  test_bdd
  ), class="berlin_data_list")

test_bddi <- structure(
  list(
    description = 'foo',
    title = 'bar',
    link = './data/data-datasetpage2.html'
  ), class = "berlin_data_dataset_info")

test_bl2 <- structure(list(
  test_bddi,
  test_bddi
), class="berlin_data_list")

# 2 different kinds of berlin_data_list:
# test_bl1 contains a list of objects with class berlin_data_dataset
# test_bl2 contains a list of objects with class berlin_data_dataset_info


test_that("berlin-data-list gives correct output for base generic methods", {
  expect_true(is.berlin_data_list(test_bl1))
  expect_equivalent(class(test_bl1), class(test_bl2[1]))
  expect_that(class(test_bl1), not(is_equivalent_to(class(test_bl2[[1]]))))
  expect_equivalent(dim(as.data.frame(test_bl1)), c(3, 7))
  expect_equivalent(dim(as.data.frame(test_bl2)), c(2, 4))
  expect_output(summary(test_bl2), '2 datasets')
})

test_that("berlin-data-list download method gives correct output", {
  test_data <- download(test_bl1)
  expect_message(download(test_bl1), 'Downloading all resources for 3 datasets')
  expect_equivalent(length(test_bl1), 3)
  expect_that(download(test_bl1, message.on.succeed=FALSE), 
              not(shows_message("Downloaded XML")))
  expect_message(download(test_bl2), 'Please call getDatasetMetaData')
  expect_equal(length(download(test_bl2)), 0)
})

test_that("berlin-data-list getDatasetMetaData method gives correct output", {
  expect_equivalent(test_bl1, getDatasetMetaData(test_bl1))
  expect_that(test_bl2, not(is_equivalent_to(getDatasetMetaData(test_bl2))))
  expect_equivalent(class(getDatasetMetaData(test_bl2)[[1]]), "berlin_data_dataset")
})
