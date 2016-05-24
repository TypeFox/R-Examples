## =============================================================================
## example : the logistic growth model
## =============================================================================
require(FME)

# 1. the analytical solution
#---------------------

model <- function(parms, time)
  with (as.list(parms),
    return(data.frame(time = time, N = K/(1+(K/N0-1)*exp(-r*time)))))

# logistic growth model parameters
TT    <- seq(0, 40, by = 10)
N0    <- 0.1

# 2. the observations
#---------------------
N     <- c(0.1, 12, 93, 96, 100)
Obs<-data.frame(time = seq(0, 40, by = 10), N = N)

# show results, compared with "observations"
plot(Obs$time, Obs$N, col = "red",  pch = 16, cex = 2,
     main  =  "logistic growth", xlab = "time", ylab = "N")

# initial "guess" of parameters
parms <- c(r = 2, K = 10)

# run the model with initial guess and plot results
lines (model(parms, TT), lwd = 2, col = "green")

# 3. Model fitting
#---------------------
# cost function
ModelCost <- function(P) {
 out  <-model(P, TT)
 return(modCost(mod = out, obs = Obs))
}

# Fit in two steps: pseudo-random search finds vicinity of minimum
(Fita <- modFit(f = ModelCost, p = parms, method = "Pseudo",
                lower = c(0, 0.1), upper = c(10, 150), control = c(numiter = 500)))

# Levenberg-Marquardt locates the minimum
Fit <- modFit(f = ModelCost, p = Fita$par,
                lower = c(0, 0.1), upper = c(10, 200))

summary(Fit)

# plot best-fit model
times <- 0:40
lines(model(Fit$par, times), lwd = 2, col = "blue")
legend("right", c("initial", "fitted"), col = c("green", "blue"), lwd = 2)

# plot residuals
MC   <- ModelCost(Fit$par)
plot(MC, main = "residuals")

# 4. Markov chain Monte Carlo
#---------------------

var0 <- summary(Fit)$modVar

MCMC<- modMCMC(f = ModelCost, p = Fit$par,
                var0 = var0, wvar0 = 1, updatecov = 100)
plot(MCMC, Full = TRUE)
pairs(MCMC)


MCMC2<- modMCMC(f = ModelCost, p = Fit$par, var0 = var0, lower = c(0, 0),
                wvar0 = 1, updatecov = 100, ntrydr = 2, niter=3000)
plot(MCMC2, Full = TRUE)
pairs(MCMC2)


SR <- summary(sensRange(parInput = MCMC2$pars, func = model, time = 0:40, sensvar = "N"))
plot(SR, xlab = "time", ylab = "N", main = "MCMC ranges")
points(Obs)
