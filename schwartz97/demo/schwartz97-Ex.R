# # pkgname <- "schwartz97"
# # source(file.path(R.home("share"), "R", "examples-header.R"))
# # options(warn = 1)
# # library('schwartz97')
# # 
# # assign(".oldSearch", search(), pos = 'CheckExEnv')
# # cleanEx()
# # nameEx("schwartz97classes")
# # ### * schwartz97classes
# # 
# # flush(stderr()); flush(stdout())
# 
# ### Name: schwartz2f-class-hierarchy
# ### Title: Classes schwartz2f and schwartz2f.fit
# ### Aliases: schwartz2f-class schwartz2f.fit schwartz2f.fit-class
# ###   show.schwartz2f show.schwartz2f.fit show,schwartz2f-method
# ###   show,schwartz2f.fit-method
# ### Keywords: classes models
# 
# ### ** Examples
# 
# 
# obj <- schwartz2f() # create an object of class schwartz2f
# obj          # print it
# coef(obj)    # get coefficients 
# unclass(obj) # see the slots 
# 
# ## create an object of class schwartz2f.fit
# data(futures)
# fit.obj <- fit.schwartz2f(futures$wheat$price, futures$wheat$ttm / 260,
#                           deltat = 1 / 260, control = list(maxit = 3))
# fit.obj          # print it
# coef(fit.obj)    # get coefficients 
# unclass(fit.obj) # see the slots 
# 
# 
# 
# 
# #cleanEx()
# #nameEx("schwartz97coef-method")
# ### * schwartz97coef-method
# 
# #flush(stderr()); flush(stdout())
# 
# ### Name: coef-method
# ### Title: Extract parameters of schwartz2f objects
# ### Aliases: coef.schwartz2f coef.schwartz2f.fit coef,schwartz2f-method
# ###   coef,schwartz2f.fit-method coefficients,schwartz2f-method
# ###   coefficients,schwartz2f.fit-method
# ### Keywords: methods utilities
# 
# ### ** Examples
# 
# 
# ## coef-method for schwartz2f-objects:
# coef(schwartz2f())
# 
# ## coef-method for schwartz2f.fit-objects:
# ## Estimate parameters for soybean oil (but stop after 3 iterations).
# data(futures)
# fit.obj <- fit.schwartz2f(futures$soybean.oil$price, futures$soybean.oil$ttm / 260,
#                           deltat = 1 / 260, control = list(maxit = 3))
# coef(fit.obj)
# 
# 
# 
# 
# # cleanEx()
# # nameEx("schwartz97constructor")
# # ### * schwartz97constructor
# # 
# # flush(stderr()); flush(stdout())
# 
# ### Name: schwartz2f-constructor
# ### Title: Create schwartz2f objects
# ### Aliases: schwartz2f
# ### Keywords: models
# 
# ### ** Examples
# 
# 
# ## Initialize a 'schwartz2f' object with high convenience yield volatility:
# obj <- schwartz2f(sigmaE = 0.7)
# 
# plot(obj) # plot it
# 
# rstate(10, time = 1, s0 = obj) # generate 10 random variates.
# 
# ## Get the probability of the event 'the spot price is >= 100 and the
# ## convenience yield is >= 0':
# pstate(c(0, -Inf), c(100, 0), time = 10, s0 = obj) 
# 
# 
# 
# 
# #cleanEx()
# #nameEx("schwartz97distrfut")
# ### * schwartz97distrfut
# 
# #flush(stderr()); flush(stdout())
# 
# ### Name: distribution-futures
# ### Title: Schwartz two-factor Model: Distribution of Futures Prices
# ### Aliases: pfutures pfutures,ANY,ANY,ANY,numeric-method
# ###   pfutures,ANY,ANY,ANY,schwartz2f-method
# ###   pfutures,ANY,ANY,ANY,schwartz2f.fit-method dfutures
# ###   dfutures,ANY,ANY,ANY,numeric-method
# ###   dfutures,ANY,ANY,ANY,schwartz2f-method
# ###   dfutures,ANY,ANY,ANY,schwartz2f.fit-method qfutures
# ###   qfutures,ANY,ANY,ANY,numeric-method
# ###   qfutures,ANY,ANY,ANY,schwartz2f-method
# ###   qfutures,ANY,ANY,ANY,schwartz2f.fit-method rfutures
# ###   rfutures,ANY,ANY,ANY,numeric-method
# ###   rfutures,ANY,ANY,ANY,schwartz2f-method
# ###   rfutures,ANY,ANY,ANY,schwartz2f.fit-method
# ### Keywords: distribution models datagen
# 
# ### ** Examples
# 
# ## Create a "schwartz2f"-object
# model <- schwartz2f()
# 
# ## Probability
# pfutures(q = 10 * 3:9, time = 0.5, ttm = 2, model, lambda = 0.01)
# 
# ## Density
# dfutures(x = c(20, 40, 100), time = 0.5, ttm = 2, model, lambda = 0.01)
# 
# ## Quantile
# qfutures(p = 0.1 * 2:5, time = 0.5, ttm = 10, model, lambda = 0.01)
# 
# ## Sample
# sim <- rfutures(n = 1000, time = 0.5, ttm = 5, model, lambda = 0.01)
# 
# hist(sim, prob = TRUE)
# lines(seq(30, 300, length = 100),
#       dfutures(seq(30, 300, length = 100),
#                time = 0.5, ttm = 5, model, lambda = 0.01), col = "red")
# 
# ## At time 0 the futures price is a deterministic function of s0 and
# ## delta0. Therefore 3 times the same value is obtained:
# rfutures(3, time = 0, ttm = 1, model, lambda = 0)
# 
# 
# 
# 
# #cleanEx()
# #nameEx("schwartz97distrstate")
# ### * schwartz97distrstate
# 
# #flush(stderr()); flush(stdout())
# 
# ### Name: distribution-state
# ### Title: Schwartz two-factor Model: Distribution of the State Variables
# ### Aliases: pstate pstate,ANY,ANY,ANY,numeric-method
# ###   pstate,ANY,ANY,ANY,schwartz2f-method dstate
# ###   dstate,ANY,ANY,numeric-method dstate,ANY,ANY,schwartz2f-method qstate
# ###   qstate,ANY,ANY,numeric-method qstate,ANY,ANY,schwartz2f-method
# ### Keywords: distribution models
# 
# ### ** Examples
# 
# ## Create a "schwartz2f"-object
# model <- schwartz2f()
# 
# ## Probability
# pstate(lower = c(0, -Inf), upper = c(45, 0.01), time = 1, model)
# 
# ## Density
# dstate(x = c(50, 0.03), time = 2, model) 
# dstate(x = rbind(c(50, 0.03), c(50, 0.1)), time = 2, model) # x is a matrix 
# 
# ## Quantile
# qstate(p = 0.5, s0 = model)
# 
# ## Generate random numbers
# object <- schwartz2f(alpha = 0.05)
# samples <- rstate(1000, time = 2, object)
# ## ...and plot histograms
# par(mfrow = c(2, 1))
# hist(samples[,1])
# abline(v = mean(object, time = 2)[1], col = "red")
# hist(samples[,2])
# abline(v = mean(object, time = 2)[2], col = "red")
# 
# 
# 
# 
# # graphics::par(get("par.postscript", pos = 'CheckExEnv'))
# # cleanEx()
# # nameEx("schwartz97filter")
# # ### * schwartz97filter
# # 
# # flush(stderr()); flush(stdout())
# 
# ### Name: filter-futures
# ### Title: Schwartz two-factor Model: Filter futures data
# ### Aliases: filter.schwartz2f filter.schwartz2f,ANY,ANY,numeric-method
# ###   filter.schwartz2f,ANY,ANY,schwartz2f-method
# ###   filter.schwartz2f,ANY,ANY,schwartz2f.fit-method
# ### Keywords: iteration
# 
# ### ** Examples
# 
# 
# data(futures)
# 
# ## Estimate parameters for coffee data (stop after 20 iterations)
# fit.obj <- fit.schwartz2f(futures$coffee$price, futures$coffee$ttm / 260,
#                           deltat = 1 / 260, control = list(maxit = 20))
# 
# ## Filter the futures data to get the spot price and the convenience yield.
# filtered <- filter.schwartz2f(futures$coffee$price, futures$coffee$ttm / 260, fit.obj)
# 
# ## ...and plot it.
# par(mfrow = c(2, 1))
# plot(filtered$state[,1], ylab = "Spot price", type = "l")
# lines(futures$coffee$price[,1], col = "red") # Closest to maturity futures
# plot(filtered$state[,2], ylab = "Convenience yield", type = "l")
# abline(h = 0)
# 
# 
# 
# 
# # graphics::par(get("par.postscript", pos = 'CheckExEnv'))
# # cleanEx()
# # nameEx("schwartz97fitted")
# # ### * schwartz97fitted
# # 
# # flush(stderr()); flush(stdout())
# 
# ### Name: fitted-futures
# ### Title: Extract Model Fitted Values
# ### Aliases: fitted.schwartz2f.fit fitted,schwartz2f.fit-method
# 
# ### ** Examples
# 
# 
# data(futures)
# 
# ## Estimate parameters for lumber data. Fit only 'mu', 'sigmaS',
# ## 'sigmaE', and 'rho' (and stop after 100 iterations).
# fit.obj <- fit.schwartz2f(futures$lumber$price, futures$lumber$ttm / 260,
#                           opt.pars = c(s0 = FALSE, delta0 = FALSE, mu = TRUE, 
#                                        sigmaS = TRUE, kappa = FALSE, alpha = FALSE, 
#                                        sigmaE = TRUE, rho = TRUE, lambda = FALSE),
#                           alpha = 0, kappa = 1.2, lambda = 0,
#                           deltat = 1 / 260, control = list(maxit = 100))
# fit.obj
# 
# ## Get the fitted values
# fitted.futures <- fitted(fit.obj, futures$lumber$price, futures$lumber$ttm / 260)
# 
# par(mfrow = c(1, 3))
# ## Plot futures prices
# plot(as.ts(futures$lumber$price), plot.type = "single", ylab = "Futures prices",
#      col = gray(seq(0.1, 0.9, length = 4)))
# ## Plot fitted values
# plot(as.ts(fitted.futures), plot.type = "single", ylab = "Fitted values",
#      col = gray(seq(0.1, 0.9, length = 4)))
# ## Plot relative difference
# plot(as.ts((fitted.futures - futures$lumber$price) / fitted.futures), plot.type = "single",
#      ylab = "Relative difference", col = gray(seq(0.1, 0.9, length = 4)))
# 
# 
# 
# 
# # graphics::par(get("par.postscript", pos = 'CheckExEnv'))
# # cleanEx()
# # nameEx("schwartz97fitting")
# # ### * schwartz97fitting
# # 
# # flush(stderr()); flush(stdout())
# 
# ### Name: parameter-estimation
# ### Title: Schwartz 1997 two factor parameter estimation
# ### Aliases: fit.schwartz2f
# ### Keywords: models iteration optimize distribution
# 
# ### ** Examples
# 
# 
# data(futures)
# 
# ## Estimate parameters for wheat data.
# ## (little precision required with reltol = 1e-3)
# fit.obj <- fit.schwartz2f(futures$wheat$price, futures$wheat$ttm / 260,
#                           deltat = 1 / 260, control = list(reltol = 1e-3))
# 
# ## See how parameter values evolved during the fit
# plot(fit.obj, type = "trace.pars")
# 
# ## Not run: 
# ##D ## Do the same but with lower tolerance level:
# ##D high.precision.fit <- fit.schwartz2f(futures$wheat$price, futures$wheat$ttm / 260,
# ##D                                      control = list(maxit = 3000, reltol = 5e-8))
# ##D 
# ##D high.precision.fit
# ##D 
# ##D plot(high.precision.fit, type = "trace.pars") # ...concerning parameter evolution.
# ##D 
# ##D ## Alpha becomes nonsensically high, kappa (speed of mean-reversion
# ##D ## of the convenience yield) goes to zero. Solution: Choose different
# ##D ## initial values and/or hold kappa constant at 1.
# ##D 
# ##D constrained.fit <- fit.schwartz2f(futures$wheat$price, futures$wheat$ttm / 260,
# ##D                                   opt.pars = c(s0 = FALSE, delta0 = FALSE, mu = TRUE, 
# ##D                                                sigmaS = TRUE, kappa = FALSE, alpha = TRUE, 
# ##D                                                sigmaE = TRUE, rho = TRUE, lambda = TRUE),
# ##D                                   alpha = 0, kappa = 1, 
# ##D                                   control = list(maxit = 3000, reltol = 5e-8))
# ##D 
# ##D constrained.fit
# ##D 
# ##D plot(constrained.fit, type = "trace.pars")
# ##D 
# ##D ## The above parameters based on a fit - where kappa was held constant at 1 -
# ##D ## look more reasonable.
# ## End(Not run)
# ## These residuals should be iid standard normal
# model.resid <- resid(fit.obj, data = futures$wheat$price, ttm = futures$wheat$ttm / 260,
#                      type = "filter.std")
# acf(model.resid, na.action = na.pass)
# par(mfrow = c(3, 2))
# apply(model.resid, 2, function(x)plot(density(na.omit(x))))
# ## ...but are anything but iid standard normal.
# 
# ## ...though the fitted values look fine:
# fitted.futures <- fitted(fit.obj, futures$wheat$price, futures$wheat$ttm / 260)
# par(mfrow = c(1, 3))
# ### Plot futures prices
# plot(as.ts(futures$wheat$price), plot.type = "single", ylab = "Futures prices",
#      col = gray(seq(0.1, 0.9, length = 4)))
# ## Plot fitted values
# plot(as.ts(fitted.futures), plot.type = "single", ylab = "Fitted values",
#      col = gray(seq(0.1, 0.9, length = 4)))
# ## Plot relative difference
# plot(as.ts((fitted.futures - futures$wheat$price) / fitted.futures), plot.type = "single",
#      ylab = "Relative difference", col = gray(seq(0.1, 0.9, length = 4)))
# 
# 
# ## Try with kappa = 1, alpha = 0, and flexible standard deviationss of
# ## the measurement errors. Stop at 200 iterations.
# fit.obj.2 <- fit.schwartz2f(futures$wheat$price, futures$wheat$ttm / 260,
#                             opt.pars = c(s0 = FALSE, delta0 = FALSE, mu = TRUE, 
#                                          sigmaS = TRUE, kappa = FALSE, alpha = FALSE, 
#                                          sigmaE = TRUE, rho = TRUE, lambda = TRUE),
#                             alpha = 0, kappa = 1, opt.meas.sd = "all",
#                             deltat = 1 / 260, control = list(maxit = 200))
# 
# 
# model.resid.2 <- resid(fit.obj.2, data = futures$wheat$price, ttm = futures$wheat$ttm / 260,
#                        type = "filter.std")
# ## Residuals seem to be better:
# acf(model.resid.2, na.action = na.pass)
# par(mfrow = c(3, 2))
# apply(model.resid.2, 2, function(x)plot(density(na.omit(x))))
# 
# 
# ## The schwartz2f.fit-object 'fit.obj' can be used to do further calculations as
# ## pricing a put option on the wheat futures which matures in 1.1
# ## years. The option expires in 1 year. The strike price is the
# ## expected futures price in 1.1 years.
# priceoption("put", time = 1, Time = 1.1, K = pricefutures(1.1, fit.obj),
#             fit.obj)
# 
# 
# 
# ## Parameter estimation for weekly soybean meal data:
# ## Get Wednesday observations:
# futures.w <- rapply(futures, function(x)x[format(as.Date(rownames(x)), "%w") == 3,],
#                     classes = "matrix", how = "list")
# 
# ## Estimate soybean meal parameters (stop after 500 iterations).
# ## ttm (time-to-maturity) is divided by 260 as it is in unit of days and
# ## deltat = 1/52 because weekly price observations are used.
# ## Estimate all measurement error standard deviations (opt.meas.sd == "all"),
# ## but hold kappa, alpha, and lambda constant.
# soybean.meal.fit <- fit.schwartz2f(data = futures.w$soybean.meal$price,
#                                    ttm = futures.w$soybean.meal$ttm / 260,
#                                    opt.meas.sd = "all",
#                                    mu = 0, kappa = 1, alpha = 0.03, r = 0.04,
#                                    opt.pars = c(s0 = FALSE, delta0 = FALSE, mu = TRUE, 
#                                      sigmaS = TRUE, kappa = FALSE, alpha = FALSE, 
#                                      sigmaE = TRUE, rho = TRUE, lambda = FALSE),
#                                    deltat = 1 / 52, control = list(maxit = 500))
# 
# soybean.meal.fit 
# 
# plot(soybean.meal.fit, type = "trace.pars") # plot the parameter evolution
# 
# ## Plot real and predicted forward curves:
# par(mfrow = c(1, 2))
# futuresplot(futures.w$soybean.meal, type = "forward.curve")
# plot(soybean.meal.fit, type = "forward.curve", data = futures.w$soybean.meal$price,
#      ttm = futures.w$soybean.meal$ttm / 260)
# 
# 
# 
# 
# 
# # graphics::par(get("par.postscript", pos = 'CheckExEnv'))
# # cleanEx()
# # nameEx("schwartz97futures-data")
# # ### * schwartz97futures-data
# # 
# # flush(stderr()); flush(stdout())
# 
# ### Name: futures-data
# ### Title: Daily futures prices
# ### Aliases: futures-data futures
# ### Keywords: datasets
# 
# ### ** Examples
# 
# 
# data(futures)
# 
# ## Plot forward curves of lumber
# futuresplot(futures$lumber, type = "forward.curve")
# 
# ## Plot time to maturity of heating oil data
# futuresplot(futures$heating.oil, type = "ttm")
# 
# ## Make 'futures' weekly, take Wednesday data
# futures.w <- rapply(futures, function(x)x[format(as.Date(rownames(x)), "%w") == 3,],
#                     classes = "matrix", how = "list")
# 
# ## Make 'futures' monthly, take the 28th day of the month
# futures.m <- rapply(futures, function(x)x[format(as.Date(rownames(x)), "%d") == 28,],
#                     classes = "matrix", how = "list") 
# 
# ## Plot weekly lumber and monthly soybean data
# futuresplot(futures.w$lumber, type = "forward.curve", main = "Lumber") 
# futuresplot(futures.m$soybean, type = "forward.curve", main = "Soybean")
# 
# ## Not run: 
# ##D ## Convert to zoo-objects:
# ##D require(zoo)
# ##D futures.zoo <- rapply(futures, function(x)zoo(x, as.Date(rownames(x))),
# ##D                       classes = "matrix", how = "list")
# ##D 
# ##D ## ...and plot it nicely using plot.zoo:
# ##D plot(futures.zoo$heating.oil$ttm)
# ##D plot(futures.zoo$wheat$vol)
# ##D plot(futures.zoo$copper$oi)
# ## End(Not run)
# 
# ## Estimate soybean meal parameters (stop after 100 iterations).
# ## ttm (time-to-maturity) is divided by 260 as it is in unit of days and
# ## deltat = 1/52 because weekly price observations are used.
# soybean.meal.fit <- fit.schwartz2f(data = futures.w$soybean.meal$price,
#                                    ttm = futures.w$soybean.meal$ttm / 260,
#                                    deltat = 1 / 52, r = 0.04,
#                                    control = list(maxit = 100))
# soybean.meal.fit
# 
# 
# 
# # cleanEx()
# # nameEx("schwartz97futures-plot")
# # ### * schwartz97futures-plot
# # 
# # flush(stderr()); flush(stdout())
# 
# ### Name: futures-plot
# ### Title: Visualization of Futures Data
# ### Aliases: futuresplot
# ### Keywords: hplot
# 
# ### ** Examples
# 
# 
# data(futures)
# 
# ## Plot time to maturity of corn data
# futuresplot(futures$corn, type = "ttm")
# 
# ## Plot forward curves of wheat data since Jan 2010
# wheat.2010 <- lapply(futures$wheat,
#                      function(x)x[as.Date(rownames(x)) > "2010-01-01",])
# futuresplot(wheat.2010, type = "forward.curve")
# 
# 
# 
# 
# # cleanEx()
# # nameEx("schwartz97meanvcov")
# # ### * schwartz97meanvcov
# # 
# # flush(stderr()); flush(stdout())
# 
# ### Name: mean-vcov-methods
# ### Title: Expected value and variance-covariance
# ### Aliases: mean.schwartz2f mean-methods mean,schwartz2f-method
# ###   vcov.schwartz2f vcov-methods vcov,schwartz2f-method
# ### Keywords: utilities methods
# 
# ### ** Examples
# 
# 
#   mean(schwartz2f(mu = 0.1), time = 1)
# 
#   mean(schwartz2f(mu = 0.2), time = 0:3)
# 
#   vcov(schwartz2f(), time = 10)
# 
#   ## Plot a schwartz2f-object including means and standard deviations
#   plot(schwartz2f(sigmaE = 0.1), n = 50, time = 5, dt = 1 / 52)
# 
# 
# 
# 
# # cleanEx()
# # nameEx("schwartz97plotfit")
# # ### * schwartz97plotfit
# # 
# # flush(stderr()); flush(stdout())
# 
# ### Name: plot.fit-method
# ### Title: Plot Schwartz two-factor fit-objects
# ### Aliases: plot.schwartz2f.fit plot,schwartz2f.fit,missing-method
# ###   plot-fit-methods
# ### Keywords: hplot methods
# 
# ### ** Examples
# 
# 
# data(futures)
# 
# ## Estimate parameters for lumber data (stop after 100 iterations)
# fit.obj <- fit.schwartz2f(futures$lumber$price, futures$lumber$ttm / 260,
#                           deltat = 1 / 260,
#                           control = list(maxit = 100))
# 
# ## Plot parameter evolution
# plot(fit.obj, type = "trace.pars")
# 
# ## Plot the state variables
# plot(fit.obj, type = "state", data = futures$lumber$price,
#      ttm = futures$lumber$ttm / 260)
# 
# ## Plot fitted and real forward curves of wheat data since Jan 2010.
# lumber.1995 <- lapply(futures$lumber, function(x)x[as.Date(rownames(x)) < "2000-01-01",])
# par(mfrow = c(1, 2))
# plot(fit.obj, type = "forward.curve", data = lumber.1995$price,
#      ttm = lumber.1995$ttm / 260)
# futuresplot(lumber.1995)
# 
# ## Plot trajectories from the state variables
# plot(fit.obj, type = "sim")
# 
# 
# 
# 
# 
# # graphics::par(get("par.postscript", pos = 'CheckExEnv'))
# # cleanEx()
# # nameEx("schwartz97plotstate")
# # ### * schwartz97plotstate
# # 
# # flush(stderr()); flush(stdout())
# 
# ### Name: plot.state-method
# ### Title: Plot Schwartz two-factor trajectories
# ### Aliases: plot.schwartz2f plot,schwartz2f,missing-method
# ###   plot-state-methods
# ### Keywords: hplot methods
# 
# ### ** Examples
# 
# 
# object <- schwartz2f(s0 = 1, mu = 0.1, sigmaS = 0.2,
#                      delta0 = 0, kappa = 2, alpha = 0.05, sigmaE = 0.1,
#                      rho = 0.5)
# 
# plot(object, n = 50, time = 2, dt = 1/52)
# 
# 
# 
# 
# # cleanEx()
# # nameEx("schwartz97pricefutures")
# # ### * schwartz97pricefutures
# # 
# # flush(stderr()); flush(stdout())
# 
# ### Name: pricing-futures
# ### Title: Schwartz two-factor Model: Futures Prices
# ### Aliases: pricefutures pricefutures,ANY,numeric-method
# ###   pricefutures,ANY,schwartz2f-method
# ###   pricefutures,ANY,schwartz2f.fit-method
# ### Keywords: models derivative
# 
# ### ** Examples
# 
# 
# ## function call by atomic arguments
# forward.curve <- pricefutures(ttm = 0.2 * 1:10, s0 = 10, delta0 = 0,
#                               alpha = 0, lambda = 0.02, r = 0)
# plot(forward.curve, type = "b")
# 
# ## function call via schwartz2f-object. 
# obj <- schwartz2f(delta0 = 0, sigmaE = 1e-5) # Make convenience yield inactive
# forward.curve <- pricefutures(ttm = 0.2 * 1:10, s0 = obj, r = 0, alphaT = 0)
# plot(forward.curve, type = "b")
# 
# 
# 
# 
# 
# # cleanEx()
# # nameEx("schwartz97priceoption")
# # ### * schwartz97priceoption
# # 
# # flush(stderr()); flush(stdout())
# 
# ### Name: pricing-options
# ### Title: Schwartz two-factor Model: European Option Prices
# ### Aliases: priceoption priceoption,ANY,ANY,ANY,ANY,numeric-method
# ###   priceoption,ANY,ANY,ANY,ANY,schwartz2f-method
# ###   priceoption,ANY,ANY,ANY,ANY,schwartz2f.fit-method
# ### Keywords: models derivative
# 
# ### ** Examples
# 
# 
# ## The option expires in 0.5 years and the futures contract in 1 year.
# priceoption(type = "call", time = 0.5, Time = 1, K = 40, g0 = 50)
# 
# ## The price of a European put option on the spot which expires in 2.5
# ## years.
# priceoption(type = "put", time = 2.5, Time = 2.5, K = 900, lambda = 0.02,
#             g0 = schwartz2f(s0 = 1000))
# 
# 
# 
# 
# # cleanEx()
# # nameEx("schwartz97randstate")
# # ### * schwartz97randstate
# # 
# # flush(stderr()); flush(stdout())
# 
# ### Name: rand-state
# ### Title: Schwartz two-factor Model: Sampling from the State Variables
# ### Aliases: rstate rstate,ANY,ANY,numeric-method
# ###   rstate,ANY,ANY,schwartz2f-method simstate
# ###   simstate,ANY,ANY,numeric-method simstate,ANY,ANY,schwartz2f-method
# ### Keywords: models datagen
# 
# ### ** Examples
# 
# ## Create a "schwartz2f"-object
# model <- schwartz2f()
# 
# ## and sample from its distribution at time = 2.5.
# sim <- rstate(n = 1000, s0 = model, time = 2.5)
# par(mfrow = c(1, 2))
# hist(sim[,1], main = "Distribution of Spot Price")
# hist(sim[,2], main = "Distribution of Convenience Yield")
# 
# 
# ## Create a trajectory over a 6 years horizon sampled on a weekly basis.
# trajectory <- simstate(6 * 52, time = 6, s0 = model)
# par(mfrow = c(1, 2))
# plot(trajectory[,1], main = "Spot Price", type = "l")
# plot(trajectory[,2], main = "Convenience Yield", type = "l")
# 
# 
# 
# 
# # graphics::par(get("par.postscript", pos = 'CheckExEnv'))
# # cleanEx()
# # nameEx("schwartz97resid")
# ### * schwartz97resid
# 
# # flush(stderr()); flush(stdout())
# 
# ### Name: resid-futures
# ### Title: Extract Model Residuals
# ### Aliases: resid.schwartz2f.fit resid,schwartz2f.fit-method
# 
# ### ** Examples
# 
# 
# data(futures)
# 
# ## Estimate parameters for live.cattle data.
# ## (little precision required with reltol = 1e-3)
# fit.obj <- fit.schwartz2f(futures$live.cattle$price, futures$live.cattle$ttm / 260,
#                           deltat = 1 / 260,
#                           control = list(maxit = 100, reltol = 1e-3))
# 
# ## Standardized residuals 
# resid.std <- resid(fit.obj, data = futures$live.cattle$price, ttm =
#                    futures$live.cattle$ttm / 260, type = "filter.std")
# acf(resid.std, na.action = na.pass) # ...are not independent
# 
# ## Real differences
# resid.real <- resid(fit.obj, data = futures$live.cattle$price, ttm =
#                     futures$live.cattle$ttm / 260, type = "real")
# 
# plot(as.ts(resid.real / futures$live.cattle$price)) # ...are 'relatively' accurate.
# 
# 
# 
# # ### * <FOOTER>
# # ###
# # cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
# # grDevices::dev.off()
# # ###
# # ### Local variables: ***
# # ### mode: outline-minor ***
# # ### outline-regexp: "\\(> \\)?### [*]+" ***
# # ### End: ***
# # quit('no')
