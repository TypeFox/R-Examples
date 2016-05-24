context("inter_frailty")

set.seed(30)
data_test <- sim_frailty(30)
delta2 <- data_test$delta + 1

# one covariate
fit <- suppressWarnings(inter_frailty(data_test, data_test$L, data_test$R, data_test$delta, c("xi1"), c("xi1"), M = 10, max_n = 5, burn_in = 0))
test_that("just one covariate works on inter_frailty", {
  expect_is(fit$par,"numeric")
  expect_is(fit$stop_c,"numeric")
})

# different covariate for each predictor
fit <- suppressWarnings(inter_frailty(data_test, data_test$L, data_test$R, data_test$delta, c("xi1"), c("xi2"), M = 10, max_n = 5, burn_in = 0))
test_that("program works for different covariates for each predictor", {
  expect_is(fit$par,"numeric")
  expect_is(fit$stop_c,"numeric")
})

# parallelism
# cl <- snow::makeCluster(2,type="SOCK")
# doSNOW::registerDoSNOW(cl)
# fit2 <- suppressWarnings(inter_frailty(data_test, data_test$L, data_test$R, data_test$delta, c("xi1"), c("xi1"),
#                      M = 10, max_n = 4, burn_in = 1, par_cl = cl))
#
# test_that("parallel mode is working", {
#   expect_is(fit2$par,"numeric")
#   expect_is(fit2$stop_c,"numeric")
# })
#
# snow::stopCluster.default(cl)

# M = 1 gives error
test_that("M = 1 gives error", {
  expect_error(suppressWarnings(inter_frailty(data_test, data_test$L, data_test$R, data_test$delta, c("xi1"), c("xi1"), M = 1, max_n = 5, burn_in = 0)))
})

# M = 250 works
# test_that("works for M = 250", {
#   expect_is(suppressWarnings(inter_frailty(data_test, data_test$L, data_test$R, data_test$delta, c("xi1"), c("xi1"), M = 250, max_n = 5, burn_in = 0))$par, "numeric")
# })

# stops on non 0 or 1 delta
test_that("throws error if delta is not binary", {
  expect_error(inter_frailty(data_test, data_test$L, data_test$R, delta2, c("xi1"), c("xi1"), M = 10, max_n = 5, burn_in = 0))
})

# stops when covariate names are wrong
test_that("stops when covariate names are wrong", {
  expect_error(suppressWarnings(inter_frailty(data_test, data_test$L, data_test$R, data_test$delta, c("ok"), c("xi1"), M = 1, max_n = 5, burn_in = 0)))
  expect_error(suppressWarnings(inter_frailty(data_test, data_test$L, data_test$R, data_test$delta, c("xi1"), c("wrong"), M = 1, max_n = 5, burn_in = 0)))
})

# stops when negative times are found
test_that("stops when negative times are found", {
  expect_error(suppressWarnings(inter_frailty(data_test, -data_test$L, data_test$R, data_test$delta, c("ok"), c("xi1"), M = 1, max_n = 5, burn_in = 0)))
  L_neg <- data_test$L
  L_neg[1] <- -5
  expect_error(suppressWarnings(inter_frailty(data_test, L_neg, data_test$R, data_test$delta, c("xi1"), c("wrong"), M = 1, max_n = 5, burn_in = 0)))
  expect_error(suppressWarnings(inter_frailty(data_test, data_test$L, -data_test$R, data_test$delta, c("xi1"), c("wrong"), M = 1, max_n = 5, burn_in = 0)))
})

# stops if L > R
test_that("stops if L > R", {
  expect_error(suppressWarnings(inter_frailty(data_test, data_test$R, data_test$L, data_test$delta, c("ok"), c("xi1"), M = 1, max_n = 5, burn_in = 0)))
})
