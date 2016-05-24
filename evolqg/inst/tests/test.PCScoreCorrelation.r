test_that("PCScoreScorrelation returns resonable results",
{
  means <- array(rnorm(40*10), c(10, 40)) 
  cov.matrix <- RandomMatrix(40, 1, 1, 10)
  taxons <- letters[1:10]
  means.list <- alply(means, 1)
  names(means.list) <- taxons
  test.array <- PCScoreCorrelation(means, cov.matrix, taxons, FALSE)
  test.list <- PCScoreCorrelation(means.list, cov.matrix, taxons, FALSE)
  expect_equal(test.array, test.list)
  test.array.plots <- PCScoreCorrelation(means, cov.matrix, taxons, TRUE)
  test.list.plots <- PCScoreCorrelation(means.list, cov.matrix, taxons, TRUE)
  #expect_equal(test.array.plots, test.list.plots)
  expect_is(test.array.plots, "list")
  expect_is(test.array.plots[[1]], "matrix")
  expect_is(test.array.plots[[2]], "list")
  expect_equal(length(test.array.plots[[2]]), 36)
  expect_equal(names(test.array.plots), c("correlation_p.value", "plots"))
})