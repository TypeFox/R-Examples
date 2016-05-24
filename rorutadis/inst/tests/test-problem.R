context("Building problem")

test_that("problem is built correctly", {
  perf = matrix(c(200, 100, 400, 200, 100, 500, 400, 400,
                  0.2, 0.3, 0.1, 0.4, 0.6, 0.4, 0.5, 0.6), nrow = 8)
  nrClasses = 4
  strictVF = FALSE
  characteristicPoints = c(2, 3)
  criteria = c('g', 'g')
  problem <- buildProblem (perf, nrClasses, strictVF, criteria, characteristicPoints)
  
  expect_that(problem$perf, equals(perf))
  expect_that(problem$nrClasses, equals(nrClasses))
  expect_that(problem$strictVF, equals(strictVF))
  expect_that(problem$criteria, equals(criteria))
  expect_that(problem$characteristicPoints, equals(characteristicPoints))
  
  expect_true(is.null(problem$assignmentsLB))
  expect_true(is.null(problem$assignmentsUB))
  expect_true(is.null(problem$pairCompMinAssignments))
  expect_true(is.null(problem$pairCompMaxAssignments))
  expect_true(is.null(problem$desiredClassCardinalitiesLB))
  expect_true(is.null(problem$desiredClassCardinalitiesUB))
  
  strictVF = TRUE
  problem <- buildProblem (perf, nrClasses, strictVF, criteria, characteristicPoints)
  expect_that(problem$strictVF, equals(strictVF))
  
  expect_error(buildProblem (perf, 1, strictVF, criteria, characteristicPoints))
})

test_that("cannot build problem with 1 class", {
  perf = matrix(c(200, 100, 400, 200, 100, 500, 400, 400,
                  0.2, 0.3, 0.1, 0.4, 0.6, 0.4, 0.5, 0.6), nrow = 8)
  nrClasses = 4
  strictVF = FALSE
  characteristicPoints = c(1, 2)
  criteria = c('g', 'g')
  
  expect_error(buildProblem (perf, 1, strictVF, criteria, characteristicPoints))
})

test_that("length of characteristic Points equals number of criteria", {
  perf = matrix(c(200, 100, 400, 200, 100, 500, 400, 400,
                  0.2, 0.3, 0.1, 0.4, 0.6, 0.4, 0.5, 0.6), nrow=8)
  nrClasses = 4
  strictVF = FALSE
  characteristicPoints_tooShort = c(2)
  characteristicPoints_tooLong = c(2,2,3)
  criteria = c('g', 'g')
  
  expect_error(buildProblem (perf, nrClasses, strictVF, criteria, characteristicPoints_tooShort))
  expect_error(buildProblem (perf, nrClasses, strictVF, criteria, characteristicPoints_tooLong))
})

context("Editing restrictions")

