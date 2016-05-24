###################################################################
# ----------------------------------------------------------------#
# Tests for the generic function getExpectation                   #
# ----------------------------------------------------------------#
###################################################################
context("getExpectation")

# getExpectation
test_that("tests of equality of getExpectation", {
  # first scenario only normEndpoint object (K=2)
  randSeq <- genSeq(rpbrPar(rb = 2, N = 12))
  randSeq@M <- matrix(rep(c(0, 1), 6), nrow = 1)
  endp <- normEndp(mu = c(0, 1), sigma = c(1, 1))
  expect_equal(getExpectation(randSeq, endp = endp), 
                   matrix(rep(c(0, 1), 6), nrow = 1))
  # second scenario normEndpoint and selection bias
  biasSB <- selBias("CS", 2, "exact")
  expect_equal(getExpectation(randSeq, biasSB, endp), 
               matrix(rep(c(0, -1), 6), nrow = 1))
  # third scenario only selection bias
  expect_equal(getExpectation(randSeq, biasSB), 
               matrix(rep(c(0, -2), 6), nrow = 1))
  # fourth scenario normEndPoint and chronological bias
  biasCB <- chronBias("linT", 1, "exact")
  expect_equal(getExpectation(randSeq, biasCB, endp), 
               matrix(rep(c(0, 1), 6) + 1:12, nrow = 1))
  # fifth scenario only chronological bias
  expect_equal(getExpectation(randSeq, biasCB), 
               matrix(1:12, nrow = 1))
  # sixth scenario normEndpoint object (K=4)
  randSeqK4 <- genSeq(rpbrPar(rb = 4, N = 12, K = 4))
  randSeqK4@M <- matrix(rep(0:3, 3), nrow = 1)
  endpK4 <- normEndp(mu = 1:4, sigma = rep(1, 4))
  expect_equal(getExpectation(randSeqK4, endp = endpK4), 
               matrix(rep(1:4, 3), nrow = 1))
  }
)
