
context("Simulator object creation")

test_that("Argument checking works", {
  with_mock(requireNamespace = function(package, ..., quietly=FALSE) FALSE,
            expect_error(create_simulator(),
                         regexp="The pomp package is needed"))
  foo <- create_simulator()
  params <- c(gamma=24, mu=0.014, d=0.014, eta=1e-4, beta=24e-2,
              rho=0.9, S_0=-1, I_0=0, R_0=0, N_0=1e2)
  expect_error(pomp::simulate(foo, params=params))
  params <- c(gamma=24, mu=0.014, d=0.014, eta=1e-4, beta=24e-2,
              rho=-9, S_0=1, I_0=0, R_0=0, N_0=1e2)
  expect_error(pomp::simulate(foo, params=params))
})

context("Gillespie direct method simulator")

test_that(paste("Mean and stddev of stationary model over time",
                "consistent with ensemble mean and stdev of dizzy",
                "progam's implementation"), {

  params <- c(gamma=24, mu=0.014, d=0.014, eta=1e-4, beta=24e-2,
              rho=0.9, S_0=1, I_0=0, R_0=0, N_0=1e2)
  covar <- data.frame(gamma_t=c(0, 0), mu_t=c(0, 0), d_t=c(0, 0), eta_t=c(0, 0),
                      beta_t=c(0, 0), time=c(0, 1e6))
  times <- seq(0, 1e6, by=1)

  sim <- create_simulator(params=params, times=times, covar=covar,
                          process_model="SIS")
  so <- pomp::simulate(sim, as.data.frame=TRUE, seed=200)
  expect_lt(abs(mean(so[, "I"]) - 0.014375), 0.5)
  expect_lt(abs(mean(so[, "S"]) - 99.9888), 1)
  expect_lt(abs(sd(so[, "I"]) - 0.5130), 0.25)
  expect_lt(abs(sd(so[, "S"]) - 9.9963), 0.1)
})

test_that(paste("Means and final stddev of time-dependent model",
                "consistent with ensemble mean and stdev of dizzy",
                "progam's implementation"), {

  params <- c(gamma=24, mu=0.014, d=0.014, eta=1e-4, beta=0e-2,
              rho=0.9, S_0=1, I_0=0, R_0=0, N_0=1e2)
  covar <- data.frame(gamma_t=c(0, 0), mu_t=c(0, 0), d_t=c(0, 0), eta_t=c(0, 0),
                      beta_t=c(0, 3 * 24e-2), time=c(0, 60))
  times <- seq(0, 50, len=100)

  sim <- create_simulator(params=params, times=times, covar=covar,
                          process_model="SIS")
  so <- pomp::simulate(sim, as.data.frame=TRUE, seed=200, nsim=1000)
  ens_infected <- unstack(so, I~sim)
  ens_susceptible <- unstack(so, S~sim)
  dzout <- read.csv(file.path("dizzy", "out-linear-trend.csv"), nrows=100)
  expect_lt(sqrt(mean((dzout$I - rowMeans(ens_infected)) ^ 2)), 1)
  expect_lt(sqrt(mean((dzout$S - rowMeans(ens_susceptible)) ^ 2)), 1)
  dzout_fluc <- read.csv(file.path("dizzy", "out-linear-trend.csv"), skip=101,
                         header=FALSE)
  tfsds <- dzout_fluc[, 2]
  names(tfsds) <- dzout_fluc[, 1]
  expect_equal(sd(ens_infected[100, ]), tfsds["I"],
               check.attributes=FALSE, tol=0.1)
  expect_equal(sd(ens_susceptible[100, ]), tfsds["S"],
               check.attributes=FALSE, tol=0.1)

  params <- c(gamma=24, mu=0.014, d=0.014, eta=0.1, beta=0e-2,
              rho=0.9, S_0=1, I_0=0, R_0=0, N_0=1e2)
  times <- seq(0, 20, len=100)
  tf <- seq(0, 20, len=1e3)

  covar <- list()
  covar$eta_t <- params["eta"] * (1 + sin(tf / 10))
  covar$beta_t <- 12 * tf / 100
  covar$gamma_t <- params["gamma"] * sin(tf / 2)
  covar$mu_t <- tf / 100
  covar$d_t <- params["d"] * exp(-tf)
  covar$time <- tf
  covar <- as.data.frame(covar)
  sim <- create_simulator(params=params, times=times, covar=covar,
                          process_model="SIS")
  so <- pomp::simulate(sim, as.data.frame=TRUE, seed=200, nsim=100)
  ens_infected <- unstack(so, I~sim)
  ens_susceptible <- unstack(so, S~sim)
  dzout <- read.csv(file.path("dizzy", "out-multiple-moving-parameters.csv"),
                    nrows=100)
  expect_lt(sqrt(mean((dzout$I - rowMeans(ens_infected)) ^ 2)), 3)
  expect_lt(sqrt(mean((dzout$S - rowMeans(ens_susceptible)) ^ 2)), 2)
  dzout_fluc <- read.csv(file.path("dizzy",
                                   "out-multiple-moving-parameters.csv"),
                         skip=101, header=FALSE)
  tfsds <- dzout_fluc[, 2]
  names(tfsds) <- dzout_fluc[, 1]
  expect_equal(sd(ens_infected[100, ]), tfsds["I"],
               check.attributes=FALSE, tol=0.25)
  expect_equal(sd(ens_susceptible[100, ]), tfsds["S"],
               check.attributes=FALSE, tol=0.25)
})

test_that(paste("Fluctuations for large system sizes approximate AR process",
                "given by linear noise approximation"), {

  params <- c(gamma=16.59091, mu=0.02, d=0.02, eta=2e-4, beta=100e-6,
              rho=0.9, S_0=0.165779, I_0=0.001004, R_0=0.833216, N_0=1e6)
  times <- seq(0, 1000, len=1000)
  sim <- create_simulator(params=params, times=times, process_model="SIR")
  so <- pomp::simulate(sim, as.data.frame=TRUE, seed=202)
  sts <- so[, c("I", "S")]
  sim_sigma <- cov(sts)
  expect_equal(sim_sigma["I", "I"], 0.1095306 * params["N_0"],
               check.attributes=FALSE, tolerance=0.2)
  expect_equal(sim_sigma["S", "S"], 18.0094688 * params["N_0"],
               check.attributes=FALSE, tolerance=0.2)
  expect_equal(sim_sigma["I", "S"], -0.1298541 * params["N_0"],
               check.attributes=FALSE, tolerance=0.2)

  sim_ac1_II <- cov(sts[-1, "I"], sts[-nrow(sts), "I"]) / sim_sigma["I", "I"]
  sim_ac1_IS <- (cov(sts[-1, "I"], sts[-nrow(sts), "S"])
                   / prod(sqrt(diag(sim_sigma))))
  sim_ac1_SI <- (cov(sts[-1, "S"], sts[-nrow(sts), "I"])
                   / prod(sqrt(diag(sim_sigma))))
  sim_ac1_SS <- cov(sts[-1, "S"], sts[-nrow(sts), "S"]) / sim_sigma["S", "S"]
  expect_equal(sim_ac1_II, 0.2037399, tol=0.2)
  expect_equal(sim_ac1_IS, 0.8501284, tol=0.2)
  expect_equal(sim_ac1_SI,-0.9121943, tol=0.2)
  expect_equal(sim_ac1_SS, 0.3079941, tol=0.2)
})
