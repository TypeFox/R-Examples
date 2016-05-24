
## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(out.extra='style="display:block; margin: auto"', fig.align="center", tidy=FALSE)


## ----ricker_constr, results='hide'---------------------------------------
library(synlik)
ricker_sl <- synlik(simulator = rickerSimul,
                    summaries = rickerStats,
                    param = c( logR = 3.8, logSigma = log(0.3), logPhi = log(10) ),
                    extraArgs = list("nObs" = 50, "nBurn" = 50)
)


## ----ricker_simul--------------------------------------------------------
ricker_sl@data <- simulate(ricker_sl, nsim = 1, seed = 54)


## ----ricker_plot---------------------------------------------------------
ricker_sl@plotFun <- function(input, ...) plot(drop(input), type = 'l', ylab = "Pop", xlab = "Time", ...)
plot(ricker_sl)


## ----ricker_simul_stats--------------------------------------------------
tmp <- simulate(ricker_sl, nsim = 10)
dim(tmp)


## ------------------------------------------------------------------------
ricker_sl@extraArgs$obsData <- ricker_sl@data


## ------------------------------------------------------------------------
tmp <- simulate(ricker_sl, nsim = 2, stats = TRUE)
tmp


## ----, results='hide'----------------------------------------------------
checkNorm(ricker_sl)


## ----ricker_slik---------------------------------------------------------
slik(ricker_sl, 
     param  = c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)),
     nsim   = 1e3)


## ----ricker_slice--------------------------------------------------------
slice(object = ricker_sl, 
      ranges = list("logR" = seq(3.5, 3.9, by = 0.01),
                    "logPhi" = seq(2, 2.6, by = 0.01),
                    "logSigma" = seq(-2, -0.5, by = 0.02)), 
      param = c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)), 
      nsim = 1000)


## ----ricker_slice_2D-----------------------------------------------------
slice(object = ricker_sl, 
      ranges = list("logR" = seq(3.5, 3.9, by = 0.02),
                    "logPhi" = seq(2, 2.6, by = 0.02)), 
      pairs = TRUE,
      param = c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)), 
      nsim = 1000, 
      multicore = TRUE,
      ncores = 2)


## ----ricker_smcmc--------------------------------------------------------
ricker_sl <- smcmc(ricker_sl, 
                   initPar = c(3.2, -1, 2.6),
                   niter = 10, 
                   burn = 3,
                   priorFun = function(input, ...) sum(input), 
                   propCov = diag(c(0.1, 0.1, 0.1))^2, 
                   nsim = 500)


## ----ricker_continue-----------------------------------------------------
ricker_sl <- continue(ricker_sl, niter = 10)


## ----ricker_plot_smcmc---------------------------------------------------
data(ricker_smcmc)
addline1 <- function(parNam, ...) abline(h = ricker_smcmc@param[parNam], lwd = 2, lty = 2, col = 3) 
addline2 <- function(parNam, ...) abline(v = ricker_smcmc@param[parNam], lwd = 2, lty = 2, col = 3)

plot(ricker_smcmc, addPlot1 = "addline1", addPlot2 = "addline2")


## ----blow_constr---------------------------------------------------------
blow_sl <- synlik(simulator = blowSimul,
                  summaries = blowStats,
                  param = log( c( "delta" = 0.16, "P" = 6.5, "N0" = 400, 
                                  "var.p" = 0.1, "tau" = 14, "var.d" = 0.1)  ),
                  extraArgs = list("nObs" = 200, "nBurn" = 200, "steps" = 1),
                  plotFun = function(input, ...){ 
                              plot(drop(input), type = 'l', ylab = "Pop", xlab = "Time", ...)
                           }
)


## ----blow_simul----------------------------------------------------------
blow_sl@data <- simulate(blow_sl, seed = 84)
blow_sl@extraArgs$obsData <- blow_sl@data


## ----blow_smcmc----------------------------------------------------------
blow_sl <- smcmc(blow_sl, 
                 initPar = log( c( "delta" = 0.1, "P" = 8, "N0" = 600, 
                                   "sig.p" = 0.2, "tau" = 17, "sig.d" = 0.3)  ),
                 niter = 2, 
                 burn = 0,
                 propCov = diag(rep(0.001, 6)),
                 nsim = 500, 
                 prior = function(input, ...){
                           sum(input) +
                           dunif(input[4], log(0.01), log(1), log = TRUE) +
                           dunif(input[6], log(0.01), log(1), log = TRUE)
                 },
                 targetRate = 0.15,
                 multicore = FALSE
)


## ----blow_plot-----------------------------------------------------------
data(blow_smcmc)
tmpTrans <- rep("exp", 6)
names(tmpTrans) <- names(blow_smcmc@param)
plot(blow_smcmc, trans = tmpTrans)


## ----bf1-----------------------------------------------------------------
data(bf1)
blow_sl@data <- bf1$pop
blow_sl@extraArgs$obsData <- blow_sl@data


## ----stableSimul---------------------------------------------------------
stableSimul <- function(param, nsim, extraArgs, ...)
{
  if( !is.loaded("stabledist") ) library(stabledist) 
  
  # Some sanity check
  if( !( c("nObs") %in% names(extraArgs) ) ) stop("extraArgs should contain nObs")
  nObs <- extraArgs$nObs 
  stopifnot( length(param) == 4 )
  param[c(1, 3)] <- exp(param[c(1, 3)])
  if(abs(param[1] - 1) < 0.01) stop("alpha == 1 is not allowed")
    
  # Actual simulation
  output <- rstable(nObs * nObs, alpha = param[1], beta = param[2], gamma = param[3], delta = param[4])
  
  return( matrix(output, nsim, nObs) )
}


## ----stableStats---------------------------------------------------------
stableStats <- function(x, extraArgs, ...){  
  
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  
  X0 <- t( apply(x, 1, quantile, probs = seq(0.1, 0.9, length.out = 10)) )
  
  unname(X0)
}


## ----stable_constr-------------------------------------------------------
stable_sl <- synlik( simulator = stableSimul,
                     summaries = stableStats,
                     param = c(alpha = log(1.5), beta = 0.1, gamma = log(1), delta = 2),
                     extraArgs = list("nObs" = 1000),
                     plotFun = function(input, ...) hist(input, xlab = "x", main = "The data", ...)
                     )
stable_sl@data <- simulate(stable_sl, seed = 67)
plot(stable_sl)


## ----stable_smcmc--------------------------------------------------------
stable_sl <- smcmc(stable_sl, 
                   initPar = c(alpha = log(1.7), beta = -0.1, gamma = log(1.3), delta = 1.5),
                   niter = 2,
                   burn = 0,
                   priorFun = function(input, ...) { 
                                 dunif(input[1], log(1), log(2), log = TRUE) + 
                                 dunif(input[2], -1, 1, log = TRUE) 
                                 }, 
                   propCov = diag(c(0.1, 0.1, 0.1, 0.1))^2, 
                   targetRate = 0.25,
                   nsim = 200)
# plot(stable_sl, trans = c("alpha" = "exp", "gamma" = "exp"))


## ------------------------------------------------------------------------
slice(object = stable_sl, 
      ranges = list("alpha" = log(seq(1.2, 1.9, by = 0.05)), 
                    "beta"  = seq(-0.5, 0.5, by = 0.05),
                    "gamma" = log(seq(0.5, 1.9, by = 0.05)),
                    "alpha" = seq(1.1, 2.2, by = 0.05)), 
      param = stable_sl@param, 
      trans = list("alpha" = "exp", "gamma" = "exp"),
      nsim = 1000, 
      multicore = TRUE,
      ncores = 2)


