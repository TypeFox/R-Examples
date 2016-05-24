context("inter_frailty_cl")

set.seed(30)
data_test_cl <- sim_frailty_cl(100, nclus = 2)
delta2 <- data_test_cl$delta + 1

# one covariate
fit <- suppressWarnings(inter_frailty_cl(data_test_cl, data_test_cl$L, data_test_cl$R, data_test_cl$delta, c("xi1"), c("xi1"), grp = data_test_cl$clus, M = 10, max_n = 3, burn_in = 0))
test_that("just one covariate works on inter_frailty_cl", {
  expect_is(fit$par,"numeric")
  expect_is(fit$stop_c,"numeric")
})

# different covariate for each predictor
fit <- suppressWarnings(inter_frailty_cl(data_test_cl, data_test_cl$L, data_test_cl$R, data_test_cl$delta, c("xi1"), c("xi2"), grp = data_test_cl$clus, M = 10, max_n = 3, burn_in = 0))
test_that("program works for different covariates for each predictor", {
  expect_is(fit$par,"numeric")
  expect_is(fit$stop_c,"numeric")
})

# parallelism
# cl <- snow::makeCluster(2,type="SOCK")
# doSNOW::registerDoSNOW(cl)
# fit2 <- suppressWarnings(inter_frailty_cl(data_test_cl, data_test_cl$L, data_test_cl$R, data_test_cl$delta, c("xi1"), c("xi1"), grp = data_test_cl$clus,
#                                        M = 10, max_n = 4, burn_in = 1, par_cl = cl))
#
# test_that("parallel mode is working", {
#   expect_is(fit2$par,"numeric")
#   expect_is(fit2$stop_c,"numeric")
# })
#
# snow::stopCluster.default(cl)

# M = 1 gives error
test_that("M = 1 gives error", {
  expect_error(suppressWarnings(inter_frailty_cl(data_test_cl, data_test_cl$L, data_test_cl$R, data_test_cl$delta, c("xi1"), c("xi1"), grp = data_test_cl$clus, M = 1, max_n = 5, burn_in = 0)))
})

# M = 250 works
# test_that("works for M = 250", {
#   expect_is(suppressWarnings(inter_frailty_cl(data_test_cl, data_test_cl$L, data_test_cl$R, data_test_cl$delta, c("xi1"), c("xi1"), grp = data_test_cl$clus, M = 250, max_n = 5, burn_in = 0))$par, "numeric")
# })

# stops on non 0 or 1 delta
test_that("throws error if delta is not binary", {
  expect_error(inter_frailty_cl(data_test_cl, data_test_cl$L, data_test_cl$R, delta2, c("xi1"), c("xi1"), grp = data_test_cl$clus, M = 10, max_n = 5, burn_in = 0))
})

# stops when covariate names are wrong
test_that("stops when covariate names are wrong", {
  expect_error(suppressWarnings(inter_frailty_cl(data_test_cl, data_test_cl$L, data_test_cl$R, data_test_cl$delta, c("ok"), c("xi1"), grp = data_test_cl$clus, M = 1, max_n = 5, burn_in = 0)))
  expect_error(suppressWarnings(inter_frailty_cl(data_test_cl, data_test_cl$L, data_test_cl$R, data_test_cl$delta, c("xi1"), c("wrong"), grp = data_test_cl$clus, M = 1, max_n = 5, burn_in = 0)))
})

# stops when negative times are found
test_that("stops when negative times are found", {
  expect_error(suppressWarnings(inter_frailty_cl(data_test_cl, -data_test_cl$L, data_test_cl$R, data_test_cl$delta, c("ok"), c("xi1"), grp = data_test_cl$clus, M = 1, max_n = 5, burn_in = 0)))
  L_neg <- data_test_cl$L
  L_neg[1] <- -5
  expect_error(suppressWarnings(inter_frailty_cl(data_test_cl, L_neg, data_test_cl$R, data_test_cl$delta, c("xi1"), c("wrong"), grp = data_test_cl$clus, M = 1, max_n = 5, burn_in = 0)))
  expect_error(suppressWarnings(inter_frailty_cl(data_test_cl, data_test_cl$L, -data_test_cl$R, data_test_cl$delta, c("xi1"), c("wrong"), grp = data_test_cl$clus, M = 1, max_n = 5, burn_in = 0)))
})

# stops if L > R
test_that("stops if L > R", {
  expect_error(suppressWarnings(inter_frailty_cl(data_test_cl, data_test_cl$R, data_test_cl$L, data_test_cl$delta, c("ok"), c("xi1"), grp = data_test_cl$clus, M = 1, max_n = 5, burn_in = 0)))
})
