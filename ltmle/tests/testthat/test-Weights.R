context("Weights")

n <- 100
W <- rnorm(n)
A1 <- rbinom(n, 1, 0.2)
Y1 <- rbinom(n, 1, plogis(W + 2*A1))
A2 <- Y2 <- rep(NA, n)
alive <- Y1 == 0
A2[alive] <- rbinom(sum(alive), 1, 0.1)
Y2[alive] <- rbinom(sum(alive), 1, plogis(W[alive] + 3*A2[alive]))
Y2[!alive] <- 1

regimes <- list(function(row) c(1, 1), function(row) c(0, 0))
data <- data.frame(W, A1, Y1, A2, Y2)
r1 <- ltmleMSM(data, Anodes=c("A1","A2"), Ynodes=c("Y1", "Y2"), final.Ynodes=c("Y1", "Y2"), survivalOutcome=TRUE, regimes=regimes, working.msm="Y~1", summary.measures=NULL, estimate.time=FALSE, msm.weights=NULL)

test_that("constant observation weights = all weights 1", {
  r.constObs <- ltmleMSM(data, Anodes=c("A1","A2"), Ynodes=c("Y1", "Y2"), final.Ynodes=c("Y1", "Y2"), survivalOutcome=TRUE, regimes=regimes, working.msm="Y~1", summary.measures=NULL, estimate.time=FALSE, observation.weights=rep(5, n), msm.weights=NULL)
  expect_equal(summary(r1)$cmat, summary(r.constObs)$cmat, tolerance=1e-4)
})

test_that("constant MSM weights = all weights 1", {
  r.constMsm1 <- ltmleMSM(data, Anodes=c("A1","A2"), Ynodes=c("Y1", "Y2"), final.Ynodes=c("Y1", "Y2"), survivalOutcome=TRUE, regimes=regimes, working.msm="Y~1", summary.measures=NULL, estimate.time=FALSE, msm.weights=matrix(123, 2, 2))
  r.constMsm2 <- ltmleMSM(data, Anodes=c("A1","A2"), Ynodes=c("Y1", "Y2"), final.Ynodes=c("Y1", "Y2"), survivalOutcome=TRUE, regimes=regimes, working.msm="Y~1", summary.measures=NULL, estimate.time=FALSE, msm.weights=matrix(123, 2, 2), observation.weights=rep(5, n))
  r.constMsm3 <- ltmleMSM(data, Anodes=c("A1","A2"), Ynodes=c("Y1", "Y2"), final.Ynodes=c("Y1", "Y2"), survivalOutcome=TRUE, regimes=regimes, working.msm="Y~1", summary.measures=NULL, estimate.time=FALSE, msm.weights=array(7, dim=c(n, 2, 2)), observation.weights=rep(5, n))
                        
  for (est in c("tmle", "iptw")) {
    expect_equal(summary(r1, est)$cmat, summary(r.constMsm1, est)$cmat, tolerance=1e-3)
    expect_equal(summary(r1, est)$cmat, summary(r.constMsm2, est)$cmat, tolerance=1e-3)
    expect_equal(summary(r1, est)$cmat, summary(r.constMsm3, est)$cmat, tolerance=1e-3)
  }
  
})

test_that("observation weights influence result", {
  obs.w <- W - min(W) + 1e-4
  r.obs <- ltmleMSM(data, Anodes=c("A1","A2"), Ynodes=c("Y1", "Y2"), final.Ynodes=c("Y1", "Y2"), survivalOutcome=TRUE, regimes=regimes, working.msm="Y~1", summary.measures=NULL, estimate.time=FALSE, observation.weights=obs.w, msm.weights=NULL)
  expect_true(plogis(r.obs$beta[1]) - plogis(r1$beta[1]) > 0.01)
})

test_that("MSM weights influence result", {
  r.msm1 <- ltmleMSM(data, Anodes=c("A1","A2"), Ynodes=c("Y1", "Y2"), final.Ynodes=c("Y1", "Y2"), survivalOutcome=TRUE, regimes=regimes, working.msm="Y~1", summary.measures=NULL, estimate.time=FALSE, msm.weights=array(W - min(W) + 1e-4, dim=c(n, 2, 2)))
  r.msm2 <- ltmleMSM(data, Anodes=c("A1","A2"), Ynodes=c("Y1", "Y2"), final.Ynodes=c("Y1", "Y2"), survivalOutcome=TRUE, regimes=regimes, working.msm="Y~1", summary.measures=NULL, estimate.time=FALSE, msm.weights=matrix(c(10, 1, 10, 1), 2, 2))
  r.msm3 <- ltmleMSM(data, Anodes=c("A1","A2"), Ynodes=c("Y1", "Y2"), final.Ynodes=c("Y1", "Y2"), survivalOutcome=TRUE, regimes=regimes, working.msm="Y~1", summary.measures=NULL, estimate.time=FALSE, msm.weights=matrix(c(1, 1, 10, 10), 2, 2))
  r.msm.empirical <- ltmleMSM(data, Anodes=c("A1","A2"), Ynodes=c("Y1", "Y2"), final.Ynodes=c("Y1", "Y2"), survivalOutcome=TRUE, regimes=regimes, working.msm="Y~1", summary.measures=NULL, estimate.time=FALSE, msm.weights="empirical")
  expect_true(plogis(r.msm1$beta[1]) - plogis(r1$beta[1]) > 0.01)
  expect_true(plogis(r.msm2$beta[1]) - plogis(r1$beta[1]) > 0.01)
  expect_true(plogis(r.msm3$beta[1]) - plogis(r1$beta[1]) > 0.01)
  expect_true(plogis(r.msm.empirical$beta[1]) - plogis(r1$beta[1]) < -0.01)
})

test_that("integer observation weights act like making copies", {
  # point estimates should be the same, but std errors are different
  
  skip_on_cran() #this seems to work ~99 times out of 100 but once in a while fails
  sampling.weight <- 4
  index <- rep(which(W > 1), each = sampling.weight - 1)
  data2 <- rbind(data, data[index, ])
  r.copies <- ltmleMSM(data2, Anodes=c("A1","A2"), Ynodes=c("Y1", "Y2"), final.Ynodes=c("Y1", "Y2"), survivalOutcome=TRUE, regimes=regimes, working.msm="Y~1", summary.measures=NULL, estimate.time=FALSE, msm.weights=NULL)
  
  observation.weights <- rep(1, n)
  observation.weights[W > 1] <- sampling.weight
  r.weights <- ltmleMSM(data, Anodes=c("A1","A2"), Ynodes=c("Y1", "Y2"), final.Ynodes=c("Y1", "Y2"), survivalOutcome=TRUE, regimes=regimes, working.msm="Y~1", summary.measures=NULL, estimate.time=FALSE, observation.weights=observation.weights, msm.weights=NULL)
  
  expect_equal(r.copies$beta, r.weights$beta, tolerance=0.001)
  expect_true(summary(r.copies)$cmat[1, "Std. Error"] / summary(r.weights)$cmat[1, "Std. Error"] < 0.95)
})