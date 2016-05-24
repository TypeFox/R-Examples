if (requireNamespace("lintr", quietly = TRUE)) {
  context("lints")
  test_that("Package Style", {
    lintr::expect_lint_free()
  })
}

set.seed(123)

context("detrending")

test_that("Mean-based detrending works", {
  expect_equal(detrend(1:10, trend="grand_mean")$x, matrix(1:10 - 5.5))
  expect_equal(detrend(cbind(1:10, 2:11), trend="grand_mean")$x,
               cbind(1:10, 2:11) - 6)
  expect_equal(detrend(1:10, trend="ensemble")$x, matrix(rep(0, 10)))
  expect_equal(detrend(cbind(1:10, 2:11), trend="ensemble")$x,
               cbind(1:10, 2:11) - 1:10 - 0.5)
})

test_that("Kernel-based detrending works", {
  x <- runif(1:10)
  x_ind <- seq_along(x)
  expect_equal(
    detrend(x, trend="local_constant", bandwidth=2, kernel="uniform")$x,
matrix(x - stats::ksmooth(x=x_ind, y=x, "box", bandwidth=2, x.points=x_ind)$y))
  expect_equal(
    detrend(x, trend="local_constant", bandwidth=3, kernel="uniform")$x,
matrix(x - stats::ksmooth(x=x_ind, y=x, "box", bandwidth=4, x.points=x_ind)$y))
})

test_that("Skipping detrending works", {
  expect_equal(detrend(1:10, trend="assume_zero")$x, matrix(1:10))
})

if (requireNamespace("np", quietly = TRUE)) {
  if (is.null(options("np.messages")$np.messages)) options(np.messages = TRUE)
  if (is.null(options("np.tree")$np.tree)) options(np.tree = FALSE)

  context("smoothing")

  test_that("Smoothing function works as expected", {
    np_smooth <- function(data, est, bandwidth, kernel="gaussian"){
      is.constant <- grepl("constant", est)
      rt <- ifelse(is.constant, "lc", "ll")
      bw <- np::npregbw(formula=rmn ~ step, bws=bandwidth,
                        regtype=rt, ckertype=kernel,
                        bandwidth.compute=FALSE, data=data,
                        na.action=na.fail)
      mod <- np::npreg(bw)
      list(smooth=fitted(mod), bandwidth=bw$bw)
    }
    data <- data.frame(step=1:10, rmn=1:10)
    expect_equal(smooth(data, est="local_constant", bandwidth=2)$smooth,
                 np_smooth(data, est="local_constant",  bandwidth=2)$smooth)
    expect_equal(smooth(data, est="local_linear", bandwidth=2)$smooth,
                 np_smooth(data, est="local_linear",  bandwidth=2)$smooth,
                 check.names=FALSE)

    expect_equal(smooth(data, est="local_constant", kernel="uniform",
                        bandwidth=3)$smooth,
                 np_smooth(data, est="local_constant", kernel="uniform",
                           bandwidth=3)$smooth)
    expect_equal(smooth(data, est="local_linear", bandwidth=3)$smooth,
                 np_smooth(data, est="local_linear", bandwidth=3)$smooth,
                 check.names=FALSE)

    data2 <- data.frame(step=20:100, rmn=rnorm(81))
    expect_equal(smooth(data2, est="local_constant", bandwidth=5)$smooth,
                 np_smooth(data2, est="local_constant", bandwidth=5)$smooth)
    expect_equal(smooth(data2, est="local_linear", bandwidth=5,
                        kernel="uniform")$smooth,
                 np_smooth(data2, est="local_linear", bandwidth=5,
                           kernel="uniform")$smooth,
                 check.names=FALSE)
    })
}

context("smoothing arguments")

test_that("invalid bandwidths lead to errors", {
  data <- data.frame(step=1:10, rmn=1:10)
  expect_error(smooth(data),
    regexp="argument \"bandwidth\" is missing, with no default")
  expect_error(smooth(data, bandwidth=c(1:20)),
    regexp="argument \"bandwidth\" must be provided as a single numeric value")
  expect_error(smooth(data, bandwidth="hmm"),
    regexp="argument \"bandwidth\" must be provided as a single numeric value")
  expect_error(smooth(data, bandwidth=0.5),
    regexp="argument \"bandwidth\" must be >= 1")
})

context("autocor")

test_that("invalid arguments lead to errors", {
  expect_error(autocor(NA))
  expect_error(autocor(letters), regexp="'x' must be numeric")
  expect_error(autocor(1:10, lag=-1), regexp="'lag' must be >= 0")
})

test_that("large bandwidth autocor estimates agree with acf", {
  w <- rnorm(1000)
  xnext <- function(xlast, w) 0.1 * xlast + w
  x <- Reduce(xnext, x=w, init=0, accumulate=TRUE)

  stats_est <- stats::acf(x, lag.max=1, plot=FALSE)$acf[2, 1, 1]
  spaero_est <- autocor(x, bandwidth=length(x) * 10)$smooth[2]
  expect_equal(stats_est, spaero_est, tolerance=0.05)

  stats_est <- stats::acf(x, lag.max=1, plot=FALSE,
                          type="covariance")$acf[2, 1, 1]
  spaero_est <- autocor(x, bandwidth=length(x) * 10,
                        cortype="covariance")$smooth[2]
  expect_equal(stats_est, spaero_est, tolerance=0.05)

  xnext <- function(xlast, w) 0.9 * xlast + w
  x <- Reduce(xnext, x=w, init=0, accumulate=TRUE)

  stats_est <- stats::acf(x, lag.max=1, plot=FALSE)$acf[2, 1, 1]
  spaero_est <- autocor(x, bandwidth=length(x) * 10)$smooth[2]
  expect_equal(stats_est, spaero_est, tolerance=0.05)
})

context("expected use of get_stats")

test_that(paste("estimate of ensemble stats consistent",
                "in case of time-dependent AR(1) model"), {
  make_time_dependent_updater <- function(f) {
    t <- 0
    updater <- function(xlast, w) {
      phi <- f(t)
      t <<- t + 1
      phi * xlast + w
    }
    return(updater)
  }
  phi_t <- function(t) {
    lambda <- min(-1  + 0.01 * t, 0)
    exp (lambda)
  }
  nensemble <- 3e3
  nobs <- 90
  x <- matrix(nrow=nobs, ncol=nensemble)
  for (i in seq(1, nensemble)){
    updater <- make_time_dependent_updater(f=phi_t)
    w <- rnorm(nobs - 1)
    init <- rnorm(n=1, sd=sqrt(1.16))
    x[, i] <- Reduce(updater, x=w, init=init, accumulate=TRUE)
  }

  est <- get_stats(x, center_trend="assume_zero", stat_bandwidth=3,
                   stat_trend="local_linear")
  lambda_ests <- log(est$stats$autocorrelation[-1])
  lambda_known <- seq(from=-1, by=0.01, len=nobs)
  var_known <-  1 / (1 - exp(2 * lambda_known))
  error_in_mean <- est$stats$mean

  expect_equal(est$centered$x + est$centered$center, x)
  expect_equal(mean(error_in_mean), 0)
  error_in_lambda <- lambda_ests - lambda_known[-1]
  expect_lt(mean(error_in_lambda ^ 2), 0.01)

  ac_error <- est$stats$autocor - exp(lambda_known)
  acov_error <- est$stats$autocov - exp(lambda_known) * var_known
  expect_lt(mean(ac_error ^ 2, na.rm=TRUE), 0.001)
  expect_lt(mean(acov_error ^ 2, na.rm=TRUE), 0.1)

  decay_time_error <- est$stats$decay_time[-1] -  -1 / lambda_known[-1]
  expect_lt((mean(decay_time_error ^ 2)), 0.5)

  expect_lt(mean( (est$stats$var - var_known) ^ 2), 0.1)
  expect_lt(mean(est$stats$skewness ^ 2), 0.01)
  expect_lt(mean( (est$stats$kurtosis - 3) ^ 2), 0.01)

  trend <- sin(2 * pi * (1:nobs) / nobs) + 2
  xx <- x + trend
  est <- get_stats(xx, center_trend="local_constant", center_bandwidth=3,
                   stat_bandwidth=3)
  lambda_ests <- log(est$stats$autocorrelation[-1])

  error_in_mean <- est$stats$mean - trend
  expect_equal(est$centered$x + est$centered$center, xx)
  expect_lt(mean(error_in_mean ^ 2), 0.01)
  error_in_lambda <- lambda_ests - lambda_known[-1]
  expect_lt(mean(error_in_lambda ^ 2), 0.01)

  ac_error <- est$stats$autocor - exp(lambda_known)
  acov_error <- est$stats$autocov - exp(lambda_known) * var_known
  expect_lt(mean(ac_error ^ 2, na.rm=TRUE), 0.001)
  expect_lt(mean(acov_error ^ 2, na.rm=TRUE), 0.1)

  decay_time_error <- est$stats$decay_time[-1] -  -1 / lambda_known[-1]
  expect_lt((mean(decay_time_error ^ 2)), 0.5)

  expect_lt(mean( (est$stats$var - var_known) ^ 2), 0.1)
  error_in_cv <- est$stats$coefficient_of_variation -  sqrt(var_known) / trend
  expect_lt(sqrt(mean( (error_in_cv) ^ 2)), 0.05)
  error_in_id <- est$stats$index_of_dispersion -  var_known / trend
  expect_lt(sqrt(mean( (error_in_id) ^ 2)), 0.2)
  expect_lt(mean(est$stats$skewness ^ 2), 0.01)
  expect_lt(mean( (est$stats$kurtosis - 3) ^ 2), 0.01)
})


test_that(paste("estimate of stats consistent",
                "in case of stationary AR(1) model"), {

  nobs <- 4e3
  w <- rnorm(nobs - 1)
  phi <- 0.1
  xnext <- function(xlast, w) phi * xlast + w
  x <- Reduce(xnext, x=w, init=0, accumulate=TRUE)

  est <- get_stats(x, center_trend="assume_zero", stat_bandwidth=nobs,
                   stat_kernel="uniform", stat_trend="local_linear")
  lambda_ests <- log(est$stats$autocorrelation[-1])
  lambda_known <- rep(log(phi), nobs)
  var_known <-  1 / (1 - exp(2 * lambda_known))
  error_in_mean <- est$stats$mean

  expect_equal(est$centered$x[, 1] + est$centered$center, x)
  expect_equal(mean(error_in_mean), 0)
  error_in_lambda <- lambda_ests - lambda_known[-1]
  expect_lt(mean(error_in_lambda ^ 2), 0.1)

  ac_error <- est$stats$autocor - exp(lambda_known)
  acov_error <- est$stats$autocov - exp(lambda_known) * var_known
  expect_lt(mean(ac_error ^ 2, na.rm=TRUE), 0.002)
  expect_lt(mean(acov_error ^ 2, na.rm=TRUE), 0.1)

  decay_time_error <- est$stats$decay_time[-1] -  -1 / lambda_known[-1]
  expect_lt((mean(decay_time_error ^ 2)), 0.5)

  expect_lt(mean( (est$stats$var - var_known) ^ 2), 0.1)
  expect_lt(mean(est$stats$skewness ^ 2), 0.02)
  expect_lt(mean( (est$stats$kurtosis - 3) ^ 2), 0.05)

  trend <- sin(2 * pi * (1:nobs) / nobs) + 2
  xx <- x + trend
  est <- get_stats(xx, center_trend="local_constant",
                   center_bandwidth=nobs / 16,
                   stat_bandwidth=nobs)
  lambda_ests <- log(est$stats$autocorrelation[-1])

  error_in_mean <- est$stats$mean - trend
  expect_equal(est$centered$x[, 1] + est$centered$center, xx)
  expect_lt(mean(error_in_mean ^ 2), 0.05)
  error_in_lambda <- lambda_ests - lambda_known[-1]
  expect_lt(mean(error_in_lambda ^ 2), 0.2)

  ac_error <- est$stats$autocor - exp(lambda_known)
  acov_error <- est$stats$autocov - exp(lambda_known) * var_known
  expect_lt(mean(ac_error ^ 2, na.rm=TRUE), 0.002)
  expect_lt(mean(acov_error ^ 2, na.rm=TRUE), 0.1)

  decay_time_error <- est$stats$decay_time[-1] -  -1 / lambda_known[-1]
  expect_lt((mean(decay_time_error ^ 2)), 0.5)

  expect_lt(mean( (est$stats$var - var_known) ^ 2), 0.1)
  error_in_cv <- est$stats$coefficient_of_variation -  sqrt(var_known) / trend
  expect_lt(sqrt(mean( (error_in_cv) ^ 2)), 0.2)
  error_in_id <- est$stats$index_of_dispersion -  var_known / trend
  expect_lt(sqrt(mean( (error_in_id) ^ 2)), 0.2)
  expect_lt(mean(est$stats$skewness ^ 2), 0.01)
  expect_lt(mean( (est$stats$kurtosis - 3) ^ 2), 0.02)
})

test_that(paste("Estimate of stats consistent with other methods",
                "in case of moving window estimates in",
                "nonstationary AR(1) model"), {

  params <- c(gamma=24, mu=0.014, d=0.014, eta=1e-4, beta=0,
              rho=0.9, S_0=1, I_0=0, R_0=0, N_0=1e5)
  covar <- data.frame(gamma_t=c(0, 0), mu_t=c(0, 0), d_t=c(0, 0), eta_t=c(0, 0),
                      beta_t=c(0, 24e-5), time=c(0, 300))
  times <- seq(0, 200, by=1 / 12)

  sim <- create_simulator(params=params, times=times, covar=covar)
  so <- pomp::simulate(sim, as.data.frame=TRUE, seed=272)

  bw <- 720
  n <- nrow(so)
  sp <- get_stats(diff(so[, "reports"]), center_kernel="uniform",
                  center_trend="local_constant", center_bandwidth=bw,
                  stat_bandwidth=bw, stat_kernel="uniform")
  ew <- earlywarnings::generic_ews(diff(so[, "reports"]),
                                   winsize=2 * bw / n * 100, detrending="no")
  spm <- lapply(sp$stats, function(x) x[(bw):(n - bw)])

  expect_equal(ew$acf1, spm$autocorrelation, tolerance=0.01)
  expect_equal(ew$sd, sqrt(spm$variance), tolerance=0.01)
  expect_equal(ew$kurt, spm$kurtosis, tolerance=0.01)
  expect_equal(ew$sk, abs(spm$skewness), tolerance=0.06)
})
