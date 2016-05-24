## ------------------------------------------------------------------------
library(spaero)
sim <- create_simulator()
simout <- pomp::simulate(sim)

## ------------------------------------------------------------------------
as(simout, "data.frame")

## ------------------------------------------------------------------------
pars <- sim@params
pars["rho"] <- 0.5
pomp::simulate(sim, params=pars, times=seq(1, 4), as.data.frame=TRUE)

## ------------------------------------------------------------------------
pomp::simulate(sim, nsim=2, times=seq(1, 2), as.data.frame=TRUE)

## ------------------------------------------------------------------------
pomp::simulate(sim, seed=342, as.data.frame=TRUE)
pomp::simulate(sim, seed=342, as.data.frame=TRUE)

## ------------------------------------------------------------------------
sim_sis <- create_simulator(process_model="SIS")
pomp::simulate(sim_sis, as.data.frame=TRUE)

## ------------------------------------------------------------------------

params <- c(gamma=24, mu=0.014, d=0.014, eta=1e-4, beta=0,
            rho=0.9, S_0=1, I_0=0, R_0=0, N_0=1e5)
covar <- data.frame(gamma_t=c(0, 0), mu_t=c(0, 0), d_t=c(0, 0), eta_t=c(0, 0),
                    beta_t=c(0, 24e-5), time=c(0, 300))
times <- seq(0, 200, by=1/12)

sim <- create_simulator(params=params, times=times, covar=covar)
so <- pomp::simulate(sim, as.data.frame=TRUE, seed=272)
plot(ts(so[, "reports"], freq=12), ylab="No. reports")


## ------------------------------------------------------------------------

st1 <- get_stats(so[, "reports"], center_kernel="uniform",
                 center_trend="local_constant", center_bandwidth=360,
                 stat_bandwidth=360)
plot_st <- function(st) {
  plot_vars <- ts(cbind(Residual=st$centered$x[, 1], Mean=st$stats$mean,
                  Autocorrelation=st$stats$autocor, Variance=st$stats$var,
                  Skewness=st$stats$skew, Kurtosis=st$stats$kurt), freq=12)
  plot(plot_vars, main="")
}
plot_st(st1)


## ------------------------------------------------------------------------

st2 <- get_stats(so[, "reports"], center_kernel="gaussian",
                 center_trend="local_constant", center_bandwidth=360,
                 stat_bandwidth=360, stat_kernel="gaussian")
plot_st(st2)



## ------------------------------------------------------------------------

st3 <- get_stats(so[, "reports"], center_kernel="gaussian",
                 center_trend="local_constant", center_bandwidth=720,
                 stat_bandwidth=720, stat_kernel="gaussian")
plot_st(st3)


