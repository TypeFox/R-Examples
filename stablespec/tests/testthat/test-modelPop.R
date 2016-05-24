context("Expected input arguments and output in modelPop")

test_that("Incorrect/missing input arguments yields errors in modelPop", {

  expect_error(modelPop(nPop="a", numVar=6, longitudinal=FALSE,
                        consMatrix=matrix(c(1, 2), 1, 2)),
               "Argument nPop should be positive numeric, e.g., 50.")

  expect_error(modelPop(nPop=25, numVar="a", longitudinal=FALSE,
                        consMatrix=matrix(c(1, 2), 1, 2)),
               "Argument numVar should be positive numeric, e.g., 6.")

  expect_error(modelPop(nPop=25, numVar=NULL, longitudinal=FALSE,
                        consMatrix=matrix(c(1, 2), 1, 2)),
               "Argument numVar cannot be missing.")

  expect_error(modelPop(nPop=25, numVar=6, longitudinal=2,
                        consMatrix=matrix(c(1, 2), 1, 2)),
               "Argument longitudinal should be either logical TRUE or FALSE.")

  expect_error(modelPop(nPop=25, numVar=6, longitudinal=NULL,
                        consMatrix=matrix(c(1, 2), 1, 2)),
               "Argument longitudinal cannot be missing.")

  expect_error(modelPop(nPop=25, numVar=6, longitudinal=FALSE,
                        consMatrix=1),
               "The constraints should be formed in a matrix.")

  expect_error(modelPop(nPop=25, numVar=6, longitudinal=FALSE,
                        consMatrix=NULL),
               "Argument consMatrix cannot be missing.")

})

test_that("Correct input arguments yield expected output in modelPop", {

  expect_true(is.matrix(modelPop(nPop=25, numVar=6, longitudinal=FALSE,
                                 consMatrix=matrix(c(1, 2), 1, 2))))

  expect_equal(nrow(modelPop(nPop=25, numVar=6, longitudinal=FALSE,
                             consMatrix=matrix(c(1, 2), 1, 2))), 25)

  expect_equal(nrow(modelPop(nPop=1, numVar=6, longitudinal=FALSE,
                             consMatrix=matrix(c(1, 2), 1, 2))), 16)

  expect_equal(nrow(modelPop(nPop=16, numVar=6, longitudinal=FALSE,
                             consMatrix=matrix(c(1, 2), 1, 2))), 16)

  expect_equal(nrow(modelPop(nPop=25, numVar=6, longitudinal=FALSE,
                             consMatrix=matrix(c(1, 2), 1, 2))), 25)

  expect_equal(nrow(modelPop(nPop=1, numVar=6, longitudinal=TRUE,
                             consMatrix=matrix(c(1, 2), 1, 2))), 52)

  expect_equal(nrow(modelPop(nPop=52, numVar=6, longitudinal=TRUE,
                             consMatrix=matrix(c(1, 2), 1, 2))), 52)

  expect_equal(nrow(modelPop(nPop=60, numVar=6, longitudinal=TRUE,
                             consMatrix=matrix(c(1, 2), 1, 2))), 60)


})
