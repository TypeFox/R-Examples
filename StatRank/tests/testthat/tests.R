
test_that("Evaluation works as expected", {
  expect_equal(Evaluation.KendallTau(c(1, 2), c(2, 1)), -1)
  expect_equal(Evaluation.KendallTau(c(1, 2), c(1, 2)), 1)
  
  expect_equal(Evaluation.LocationofWinner(c(1, 2, 3), c(1, 2, 3)), 1)
  expect_equal(Evaluation.LocationofWinner(c(1, 2, 3), c(3, 2, 1)), 3)
  
  expect_equal(TVD(c(.1, .2, .7), c(.2, .1, .7)), .1)
  expect_equal(TVD(c(.1, .2, .7), c(.1, .2, .7)), 0)
})


test_that("Helpers works as expected", {
  expect_equal(scores.to.order(c(1, 2, 3)), c(3, 2, 1))
  expect_equal(scores.to.order(c(1, 0, -.1)), c(1, 2, 3))
  
  data(Data.Test)
  Data.Test.pairs <- Breaking(Data.Test, "full")
  expect_equal(generateC(Data.Test.pairs, 5), matrix(c(0, 0.6, 0.8, 0.6, 0.8, 0.4,0, 0.6, 0.6, 0.6, 0.2, 0.4,0, 0.6, 0.6, 0.4,0.4,0.4,0, 0.6, 0.2, 0.4,0.4,0.4,0), byrow=TRUE, ncol = 5))
  expect_equal(generateC(Data.Test.pairs, 5, normalized = FALSE), matrix(c(0, 3, 4, 3, 4, 2, 0, 3, 3, 3, 1, 2, 0, 3, 3, 2, 2, 2, 0, 3, 1, 2, 2, 2, 0), byrow=TRUE, ncol = 5))
})