## Simulation study demonstrating that estimation methods work, i.e., estimate the parameters correctly
## This takes a while to run
expit <- function(x) exp(x)/(1+exp(x))

library(pseval)
library(survival)

gen_binary <- function(n){

  Z <- rbinom(n, 1, .5)
  X <- rbinom(n, 1, .35)  ## pre-treatment covariate

  ## generate S conditional on Z
  ps.0 <- expit(X + rnorm(n, sd = .1))
  ps.1 <- expit(2 + X + rnorm(n, sd = .1))

  S.0 <- rbinom(n, 1, ps.0)
  S.1 <- rbinom(n, 1, ps.1)

  risk.obs <- (1 - 0.5 * S.1 - 0 * Z - 1 * S.1 * Z)
  risk.0 <- (1 - 0.5 * S.1)
  risk.1 <- (1 - 1.5 * S.1)

  Y.0 <- rbinom(n, 1, expit(risk.0))
  Y.1 <- rbinom(n, 1, expit(risk.1))
  Y.obs <- ifelse(Z == 1, Y.1, Y.0)

  S.1[Z == 0] <- NA
  S.0[Z == 1] <- NA
  S.obs <- ifelse(Z == 1, S.1, S.0)

  data.frame(Z, BIP = X, S.obs, Y.obs)

}
## The true model for the time is exponential, with parameters (intercept) = 1, S(1) = -0.5, Z = 0, S(1):Z = -1. The true model for binary is logistic, with the same parameter values.

swoop <- function(psdesign){

  c(psdesign$estimates$par, psdesign$estimates$convergence == 0)

}

true.par <- c(1, -.5, 0, -1)
names(true.par) <- c("(Intercept)", "S.1", "Z", "S.1:Z")

nsim <- 10
est.bin.para <- matrix(NA, nrow = nsim, ncol = 5)
est.bin.semi <- matrix(NA, nrow = nsim, ncol = 5)
est.bin.pseudo <- matrix(NA, nrow = nsim, ncol = 5)
est.surv.para <- matrix(NA, nrow = nsim, ncol = 5)
est.surv.semi <- matrix(NA, nrow = nsim, ncol = 5)

if(FALSE){ # change to TRUE if you want to run simulation
for(i in 1:nsim){

  fakedata <- generate_example_data(n = 800)
  fakedat.bin <- gen_binary(n = 800)

  binary.ps <- psdesign(data = fakedata, Z = Z, Y = Y.obs, S = S.obs, BIP = BIP)
  categ.ps <- psdesign(data = fakedat.bin, Z = Z, Y = Y.obs, S = S.obs, BIP = BIP)
  surv.ps <- psdesign(data = fakedata, Z = Z, Y = Surv(time.obs, event.obs), S = S.obs, BIP = BIP, tau = 0)

  est.bin.para[i, ] <- swoop(binary.ps + integrate_parametric(S.1 ~ BIP) + risk_binary(D = 1000) + ps_estimate())
  est.bin.semi[i, ] <- swoop(binary.ps + integrate_semiparametric(S.1 ~ BIP, S.1 ~ 1) + risk_binary(D = 1000) + ps_estimate())
  est.bin.pseudo[i, ] <- swoop(categ.ps + integrate_nonparametric(S.1 ~ BIP) + risk_binary(D = 1000) + ps_estimate(method = "pseudo-score"))

  est.surv.para[i, ] <- swoop(surv.ps + integrate_parametric(S.1 ~ BIP) + risk_exponential(D = 1000) + ps_estimate())
  est.surv.semi[i, ] <- swoop(surv.ps + integrate_semiparametric(S.1 ~ BIP, S.1 ~ 1) + risk_exponential(D = 1000) + ps_estimate())

}

res <- list(bin.para = est.bin.para, bin.semi = est.bin.semi, bin.pseudo = est.bin.pseudo, surv.para = est.surv.para, surv.semi = est.surv.semi)

plotdem <- function(mat){

  for(i in 1:(ncol(mat) - 1)){
    hist(mat[, i])
    abline(v = true.par[i], col = "red", lwd = 2)
  }

  colMeans(mat[mat[, ncol(mat)] == 1, -ncol(mat)]) - true.par

}

par(mfrow = c(2, 2))
lapply(res, plotdem)

}