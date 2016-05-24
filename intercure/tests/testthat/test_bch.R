context("inter_bch")

# no cure threshold
test_that("stops if no cure threshold is found", {
  expect_error(inter_bch(test_no_cure,
                         test_no_cure$L,
                         test_no_cure$R,
                         "Trt")
)
})

# one covariate
fit <- inter_bch(data_test_frailty, data_test_frailty$L, data_test_frailty$R, "xi1")
test_that("just one covariate works on inter_bch", {
  expect_is(fit$par,"numeric")
  expect_is(fit$stop_c,"numeric")
})

# expect result (safety for changes)
test_that("results doesn't change", {
  expect_equal(fit$par, c(-0.6035598, -0.2570020), tolerance = 0.0001)
})

# stops when covariate names are wrong
test_that("stops when covariate names are wrong", {
  expect_error(inter_bch(data_test_frailty, data_test_frailty$L, data_test_frailty$R, "ok"))
  expect_error(inter_bch(data_test_frailty, data_test_frailty$L, data_test_frailty$R, c("ok","xi1")))
})

# stops when negative times are found
test_that("stops when negative times are found", {
  expect_error(inter_bch(data_test_frailty, -data_test_frailty$L, data_test_frailty$R, "xi1"))
  L_neg <- data_test_frailty$L
  L_neg[1] <- -5
  expect_error(inter_bch(data_test_frailty, L_neg, data_test_frailty$R, "xi1"))
  expect_error(inter_bch(data_test_frailty, data_test_frailty$L, -data_test_frailty$R, "xi1"))
})

# stops if L > R
test_that("stops if L > R", {
  expect_error(inter_bch(data_test_frailty, data_test_frailty$R, data_test_frailty$L, "xi1"))
})
