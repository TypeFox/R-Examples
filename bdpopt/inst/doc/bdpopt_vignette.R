### R code from vignette source 'bdpopt_vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: bdpopt_vignette.Rnw:255-275 (eval = FALSE)
###################################################
## ## Perform an optimisation for the simple normal model
## out <- n.opt(nu = 0, tau = 1, sigma = 1, alpha = 0.025,
##              gain.constant = 1, gain.function = function(X, mu) 0,
##              fixed.cost = 0, sample.cost = 0.005,
##              k = 1, n.min = 1, n.max = 50, n.step = 1, 
##              n.iter = 10000, n.burn.in = 1000, n.adapt = 1000,
##              regression.type = "loess",
##              plot.results = TRUE, independent.SE = FALSE,
##              parallel = FALSE, path.to.package = NA)   
## 
## ## Print the grid points used for the sample size,
## ## and the corresponding estimates of the expected utility
## print(out$ns)
## print(out$eus)
## 
## ## Print the estimate of the optimal sample size,
## ## and the corresponding utility
## print(out$opt.arg)
## print(out$opt.eu)
## 


###################################################
### code chunk number 2: bdpopt_vignette.Rnw:278-323 (eval = FALSE)
###################################################
## ## Mean and precision parameters for the priors
## theta.mu <- c(0, 2, 0, 0); theta.tau <- c(1, 1, 8, 8)
## eta.mu <- c(0, 2, 0, 0); eta.tau <- c(1, 1, 8, 8)    
## 
## ## Sample size and doses for each observation in phase II
## n.II <- rep(10, 10); d.II <- seq(0.1, 1, 0.1)
## 
## ## Observed responses phase II responses,
## ## taken from the efficacy and safety models using the parameter values
## ## theta = eta = c(0, 2, 1, 1) (rounded to two decimals).
## YE.II <- c(0.18, 0.33, 0.46, 0.57, 0.67, 0.75, 0.82, 0.89, 0.95, 1.00)
## YS.II <- c(0.18, 0.33, 0.46, 0.57, 0.67, 0.75, 0.82, 0.89, 0.95, 1.00)
## 
## sigmaE <- 1; sigmaS <- 1 ## Standard deviations            
## k.III <- 2 ## Number of phase III trials
##     
## m1 <- create.normal.model(theta.mu, theta.tau, eta.mu, eta.tau,
##                           n.II, d.II, YE.II, YS.II,
##                           sigmaE, sigmaS, k.III, path.to.package = NA)
## 
## ## Define a utility function
## n.min <- 0; sig.level<- 0.025; safety.max <- 0.6
## cE <- 1300; cS <- -1000; p <- 0.5
## fixed.cost <- 10; cost.per.sample <- 0.2
## 
## u <- create.utility.function(m1, n.min, sig.level, safety.max,
##                              cE, cS, p, fixed.cost, cost.per.sample)    
## 
## ## Define a grid and simulate the utility for each grid point
## n.iter <- 4000; n.burn.in <- 1000; n.adapt <- 1000
## gsl <- list(n.III = list(c(1), list(c(10, 150, 10))),
##             d.III = list(c(1), list(c(0.1, 0.4, 0.1))))    
## m2 <- eval.on.grid(m1, u, gsl, n.iter, n.burn.in, n.adapt,
##                    independent.SE = FALSE, parallel = TRUE)
##     
## ## Do gaussian process regression for the model
## m3 <- fit.gpr(m2, start = c(30, 50, 0.2), gr = TRUE, method = "L-BFGS-B",
##                   lower = c(10, 10, 0.1), upper = Inf)
## 
## ## Plot the results of the evaluation and gpr regression
## plot(m3, "n.III[1]", fixed = seq(0.1, 0.4, 0.1))        
## 
## ## Optimisation (defaulting to method "L-BFGS-B" of the optim function)
## optimise.eu(m3, start = c(100, 0.3))
## 


###################################################
### code chunk number 3: bdpopt_vignette.Rnw:326-343 (eval = FALSE)
###################################################
## ## Create a sequential decision problem object
## dp <- sequential.normal.dp(n.stages = 4, group.size = 10,
##                            tau = 1, sigma = 1,
##                            stage.cost = 0.1, final.cost = 1, final.gain = 2)
## 
## ## Solve the sequential decision problem and plot the results
## out <- optimise.sequential.normal.eu(dp = dp,
##                                      range = 8, step.size = 0.02,
##                                      prior.mean = 0,
##                                      n.sims = 1000,
##                                      plot.results = TRUE)
## 
## ## Print the optimal decision and corresponding expected utility
## ## at the first stage assuming a prior mean of 0
## print(out$opt.decision(1, 0))
## print(out$opt.utility(1, 0))
## 


