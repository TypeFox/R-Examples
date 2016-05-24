context('autovar')

testdata_raw_dataframe <- function() {
  data_matrix <- matrix(nrow = 40, ncol = 3)
  data_matrix[, ] <- runif(ncol(data_matrix) * nrow(data_matrix), 1, nrow(data_matrix))
  while (sum(is.na(data_matrix)) == 0)
    data_matrix[as.logical(round(runif(ncol(data_matrix) * nrow(data_matrix), -0.3, 0.7)))] <- NA
  colnames(data_matrix) <- c('rumination', 'happiness', 'activity')
  as.data.frame(data_matrix)
}


test_that('autovar function returns hello world', {
  with_mock(
    `parallel::clusterMap` = function(cluster, ...) {
      mapply(...)
    },
    `parallel::makeCluster` = function(...) {
      NULL
    },
    `parallel::stopCluster` = function(...) {
      NULL
    },
    expect_equal(class(autovar(testdata_raw_dataframe(),
                         selected_column_names = c('rumination',
                                                   'happiness',
                                                   'activity'),
                         imputation_iterations = 1)),
                 "list")
  )
})


test_that('nr_dummy_variables calls its subfunctions correctly', {
  varest <- list(datamat = matrix(nrow = 40, ncol = 3, dimnames = list(NULL, c('a', 'day_1', 'day_2'))))
  expect_equal(autovarCore:::nr_dummy_variables(varest), 1)
  varest <- list(datamat = matrix(nrow = 40, ncol = 3, dimnames = list(NULL, c('a', 'day_1', 'c'))))
  expect_equal(autovarCore:::nr_dummy_variables(varest), 1)
  varest <- list(datamat = matrix(nrow = 40, ncol = 3, dimnames = list(NULL, c('a', 'outlier_1', 'c'))))
  expect_equal(autovarCore:::nr_dummy_variables(varest), 1)
  varest <- list(datamat = matrix(nrow = 40, ncol = 3, dimnames = list(NULL, c('a', 'outlier_2', 'outlier_4'))))
  expect_equal(autovarCore:::nr_dummy_variables(varest), 2)
  varest <- list(datamat = matrix(nrow = 40, ncol = 3, dimnames = list(NULL, c('day_2', 'outlier_2', 'outlier_4'))))
  expect_equal(autovarCore:::nr_dummy_variables(varest), 3)
  varest <- list(datamat = matrix(nrow = 40, ncol = 3, dimnames = list(NULL, c('day_2', 'day_4', 'outlier_3'))))
  expect_equal(autovarCore:::nr_dummy_variables(varest), 2)
})

test_that('nr_dummy_variables works with 0 dummy variables', {
  varest <- list(datamat = matrix(nrow = 40, ncol = 3, dimnames = list(NULL, c('a', 'b', 'c'))))
  expect_equal(autovarCore:::nr_dummy_variables(varest), 0)
})


test_that('insert_model_into_list works when an empty list of models is given', {
  model <- list(bucket = 0.05, nr_dummy_variables = 1, model_score = 100)
  expected_result <- list(model)
  expect_equal(autovarCore:::insert_model_into_list(model, list(), TRUE),
               expected_result)
})

test_that('insert_model_into_list works when the model is to be appended to the end of the list', {
  model <- list(bucket = 0.05, nr_dummy_variables = 1, model_score = 100)
  model_list <- list(list(bucket = 0.05, nr_dummy_variables = 0, model_score = 100))
  expected_result <- append(model_list, list(model))
  expect_equal(autovarCore:::insert_model_into_list(model, model_list, TRUE),
               expected_result)
})

test_that('insert_model_into_list works when the model is to be prepended to the list', {
  model <- list(bucket = 0.05, nr_dummy_variables = 1, model_score = 100)
  model_list <- list(list(bucket = 0.05, nr_dummy_variables = 1, model_score = 101))
  expected_result <- append(model_list, list(model), after = 0)
  expect_equal(autovarCore:::insert_model_into_list(model, model_list, TRUE),
               expected_result)
})

test_that('insert_model_into_list inserts the model at the correct position', {
  model_list <- list(list(bucket = 0.05, nr_dummy_variables = 1, model_score = 103),
                     list(bucket = 0.05, nr_dummy_variables = 2, model_score = 101),
                     list(bucket = 0.01, nr_dummy_variables = 1, model_score = 99),
                     list(bucket = 0.01, nr_dummy_variables = 1, model_score = 101),
                     list(bucket = 0.01, nr_dummy_variables = 2, model_score = 101),
                     list(bucket = 0.005, nr_dummy_variables = 0, model_score = 101),
                     list(bucket = 0.005, nr_dummy_variables = 0, model_score = 102))
  model <- list(bucket = 0.01, nr_dummy_variables = 1, model_score = 100)
  expected_result <- append(model_list, list(model), after = 3)
  expect_equal(autovarCore:::insert_model_into_list(model, model_list, TRUE),
               expected_result)
  model <- list(bucket = 0.01, nr_dummy_variables = 0, model_score = 100)
  expected_result <- append(model_list, list(model), after = 2)
  expect_equal(autovarCore:::insert_model_into_list(model, model_list, TRUE),
               expected_result)
  model <- list(bucket = 0.05, nr_dummy_variables = 1, model_score = 105)
  expected_result <- append(model_list, list(model), after = 1)
  expect_equal(autovarCore:::insert_model_into_list(model, model_list, TRUE),
               expected_result)
})

test_that('insert_model_into_list handles the compare_outliers argument correctly', {
  model_list <- list(list(bucket = 0.05, nr_dummy_variables = 1, model_score = 103),
                     list(bucket = 0.05, nr_dummy_variables = 2, model_score = 101),
                     list(bucket = 0.01, nr_dummy_variables = 1, model_score = 99),
                     list(bucket = 0.01, nr_dummy_variables = 1, model_score = 101),
                     list(bucket = 0.01, nr_dummy_variables = 2, model_score = 101),
                     list(bucket = 0.005, nr_dummy_variables = 0, model_score = 101),
                     list(bucket = 0.005, nr_dummy_variables = 0, model_score = 102))
  model <- list(bucket = 0.01, nr_dummy_variables = 3, model_score = 100)
  expected_result <- append(model_list, list(model), after = 3)
  expect_equal(length(expected_result), length(model_list) + 1)
  expect_equal(autovarCore:::insert_model_into_list(model, model_list, FALSE),
               expected_result)
  model <- list(bucket = 0.01, nr_dummy_variables = 0, model_score = 100)
  expected_result <- append(model_list, list(model), after = 3)
  expect_equal(autovarCore:::insert_model_into_list(model, model_list, FALSE),
               expected_result)
  model <- list(bucket = 0.05, nr_dummy_variables = 1, model_score = 105)
  expected_result <- append(model_list, list(model), after = 2)
  expect_equal(autovarCore:::insert_model_into_list(model, model_list, FALSE),
               expected_result)
})


test_that('merge_model_lists works when one of the lists is empty', {
  model_list_a <- list(list(bucket = 0.05, nr_dummy_variables = 1, model_score = 103),
                       list(bucket = 0.05, nr_dummy_variables = 2, model_score = 101),
                       list(bucket = 0.01, nr_dummy_variables = 1, model_score = 99),
                       list(bucket = 0.01, nr_dummy_variables = 1, model_score = 101),
                       list(bucket = 0.01, nr_dummy_variables = 2, model_score = 101),
                       list(bucket = 0.005, nr_dummy_variables = 0, model_score = 101),
                       list(bucket = 0.005, nr_dummy_variables = 0, model_score = 102))
  model_list_b <- list()
  expect_equal(autovarCore:::merge_model_lists(model_list_a, model_list_b, TRUE),
               model_list_a)
  expect_equal(autovarCore:::merge_model_lists(model_list_b, model_list_a, TRUE),
               model_list_a)
})

test_that('merge_model_lists works when one of the lists should be place fully behind the other', {
  model_list_a <- list(list(bucket = 0.05, nr_dummy_variables = 1, model_score = 103),
                       list(bucket = 0.05, nr_dummy_variables = 2, model_score = 101),
                       list(bucket = 0.01, nr_dummy_variables = 1, model_score = 99),
                       list(bucket = 0.01, nr_dummy_variables = 1, model_score = 101),
                       list(bucket = 0.01, nr_dummy_variables = 2, model_score = 101),
                       list(bucket = 0.005, nr_dummy_variables = 0, model_score = 101),
                       list(bucket = 0.005, nr_dummy_variables = 0, model_score = 102))
  model_list_b <- list(list(bucket = 0.005, nr_dummy_variables = 0, model_score = 103),
                       list(bucket = 0.005, nr_dummy_variables = 1, model_score = 101),
                       list(bucket = 0.005, nr_dummy_variables = 2, model_score = 99),
                       list(bucket = 0.005, nr_dummy_variables = 2, model_score = 101),
                       list(bucket = 0.001, nr_dummy_variables = 0, model_score = 100),
                       list(bucket = 0.001, nr_dummy_variables = 0, model_score = 101),
                       list(bucket = 0.001, nr_dummy_variables = 3, model_score = 102))
  expected_result <- append(model_list_a, model_list_b)
  expect_equal(length(expected_result), length(model_list_a) + length(model_list_b))
  expect_equal(autovarCore:::merge_model_lists(model_list_a, model_list_b, TRUE),
               expected_result)
  expect_equal(autovarCore:::merge_model_lists(model_list_b, model_list_a, TRUE),
               expected_result)
})

test_that('merge_model_lists works when lists need to be interleaved', {
  model_list_a <- list(list(bucket = 0.05, nr_dummy_variables = 1, model_score = 103),
                       list(bucket = 0.01, nr_dummy_variables = 1, model_score = 99),
                       list(bucket = 0.01, nr_dummy_variables = 2, model_score = 101),
                       list(bucket = 0.005, nr_dummy_variables = 0, model_score = 102),
                       list(bucket = 0.005, nr_dummy_variables = 1, model_score = 101),
                       list(bucket = 0.005, nr_dummy_variables = 2, model_score = 101),
                       list(bucket = 0.001, nr_dummy_variables = 0, model_score = 101))
  model_list_b <- list(list(bucket = 0.05, nr_dummy_variables = 2, model_score = 101),
                       list(bucket = 0.01, nr_dummy_variables = 1, model_score = 101),
                       list(bucket = 0.005, nr_dummy_variables = 0, model_score = 101),
                       list(bucket = 0.005, nr_dummy_variables = 0, model_score = 103),
                       list(bucket = 0.005, nr_dummy_variables = 2, model_score = 99),
                       list(bucket = 0.001, nr_dummy_variables = 0, model_score = 100),
                       list(bucket = 0.001, nr_dummy_variables = 3, model_score = 102))
  expected_result <- list(list(bucket = 0.05, nr_dummy_variables = 1, model_score = 103),
                          list(bucket = 0.05, nr_dummy_variables = 2, model_score = 101),
                          list(bucket = 0.01, nr_dummy_variables = 1, model_score = 99),
                          list(bucket = 0.01, nr_dummy_variables = 1, model_score = 101),
                          list(bucket = 0.01, nr_dummy_variables = 2, model_score = 101),
                          list(bucket = 0.005, nr_dummy_variables = 0, model_score = 101),
                          list(bucket = 0.005, nr_dummy_variables = 0, model_score = 102),
                          list(bucket = 0.005, nr_dummy_variables = 0, model_score = 103),
                          list(bucket = 0.005, nr_dummy_variables = 1, model_score = 101),
                          list(bucket = 0.005, nr_dummy_variables = 2, model_score = 99),
                          list(bucket = 0.005, nr_dummy_variables = 2, model_score = 101),
                          list(bucket = 0.001, nr_dummy_variables = 0, model_score = 100),
                          list(bucket = 0.001, nr_dummy_variables = 0, model_score = 101),
                          list(bucket = 0.001, nr_dummy_variables = 3, model_score = 102))
  expect_equal(length(expected_result), length(model_list_a) + length(model_list_b))
  expect_equal(autovarCore:::merge_model_lists(model_list_a, model_list_b, TRUE),
               expected_result)
  expect_equal(autovarCore:::merge_model_lists(model_list_b, model_list_a, TRUE),
               expected_result)
})
