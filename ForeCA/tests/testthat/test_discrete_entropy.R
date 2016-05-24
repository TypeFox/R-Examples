context("discrete_entropy")
nn <- 10
# does not add up to 1
kVec <- rnorm(nn)^2 + 1
kProbs <- kVec / sum(kVec)
kUniformDistr <- rep(1 / nn, nn)

test_that("discrete_entropy computes right entropy for uniform", {
  
  # entropy of uniform is always larger than any other vector
  expect_true(discrete_entropy(kProbs) < discrete_entropy(kUniformDistr))
  # 0 times 0 is zero (remove attribute)
  expect_equal(log2(length(kUniformDistr)), c(discrete_entropy(kUniformDistr)))
  # check that base works correctly
  expect_equal(log2(length(kUniformDistr))/log2(exp(1)), 
               c(discrete_entropy(kUniformDistr, base = exp(1))))
  # 'base' attribute
  expect_true(attr(discrete_entropy(kProbs), "base") == 2)
})

test_that("entropy of certain event is 0", {
  # entropy of sure event is 0
  expect_equal(0, c(discrete_entropy(c(1, 0, 0))))
})

test_that("probabilities must be non-negative for entropy", {
  # negative values are not allowed
  expect_error(discrete_entropy(c(-0.5, 0)))
})

test_that("smooting increases entropy", {
  
  # test smoothing stuff
  expect_error(c(discrete_entropy(kProbs, prior.weight = 2)))
  expect_error(c(discrete_entropy(kProbs, prior.weight = -1)))
  
  # smoothing increases entropy
  expect_true(discrete_entropy(kProbs, prior.weight = 0.1) > 
                discrete_entropy(kProbs, prior.weight = 0))
})
