context("Expected input arguments and output in repairCyclicModel")

test_that("Incorrect/missing argument inputs yield errors in repairCyclicModel", {

  num_vars <- 6
  model_a <- round(runif(num_vars * num_vars))

  expect_error(repairCyclicModel(stringModel="a", numVar=num_vars,
                                 longitudinal=FALSE),
               "Argument numVar should be binary vector.")

  expect_error(repairCyclicModel(stringModel=NULL, numVar=num_vars,
                                 longitudinal=FALSE),
               "Argument stringModel cannot be missing.")

  expect_error(repairCyclicModel(stringModel=model_a, numVar="a",
                                 longitudinal=FALSE),
               "Argument numVar should be positive numeric, e.g., 6.")

  expect_error(repairCyclicModel(stringModel=model_a, numVar=NULL,
                                 longitudinal=FALSE),
               "Argument numVar cannot be missing.")

  expect_error(repairCyclicModel(stringModel=model_a, numVar=num_vars,
                                 longitudinal="a"),
               "Argument longitudinal should be either logical TRUE or FALSE.")


  expect_error(repairCyclicModel(stringModel=model_a, numVar=num_vars,
                                 longitudinal=NULL),
               "Argument longitudinal cannot be missing.")

})

test_that("Correct input arguments yield expected output in repairCyclicModel.", {
  num_vars <- 6
  model_a <- round(runif(num_vars * num_vars))

  expect_true(all(repairCyclicModel(stringModel=model_a, numVar=num_vars,
                                longitudinal=FALSE) %in% 0:1))

  expect_true(is.numeric(repairCyclicModel(stringModel=model_a, numVar=num_vars,
                                    longitudinal=FALSE)))

  expect_equal(length(repairCyclicModel(stringModel=model_a, numVar=num_vars,
                                        longitudinal=FALSE)), 30)


})
