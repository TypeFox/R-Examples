##
## This R code illustrates how to estimates parameters and calculate posterior model probabilities for a set of constant variance models (with different 
## lags) fitted to the full England and Wales data set.
## 
## With a few minor adjustments this code can be adapted to estimate parameters for stochastic volatility and random variance shift models, also 
## illustrated in the paper. However, the MCMC of these models will take a lot longer (especially the RV models), and hence are not demonstrated here. 
## 
## The model averaging code can also be easily adapted to handle a larger model space, once parameters and normalising constants are estimated.
##
## Code should work for alternative data sets (such as shortened series).
##

##
## Create differenced population rate from England and Wales population totals (ew in tsbugs package)
##
library("tsbugs")
r <- ts(ew[2:167]/ew[1:166] - 1, start = 1841)
y <- diff(r)

##
## Create folders for mcmc simulations for each model
## 
# this is where the folders are going...
getwd()
# create folders
for (p in 0:8) {
  dir.create(paste0("AR", p))
}

##
## Create BUGS model and estimate parameters using MCMC via R2OpenBUGS
##
library("R2OpenBUGS")
library("coda")

for (p in 0:8) {
  # BUGS model with order p
  cv <- ar.bugs(y = y, ar.order = p, beg = 9)
  # write the BUGS model to .txt file in the working directiory for use in the bugs function of R2OpenBUGS
  writeLines(cv$bug, "cv.txt")
  
  # MCMC via R2OpenBUGS
  cv.bug <- bugs(data = cv$data, 
                 inits = list(inits(cv, warn = FALSE)), 
                 param = c(nodes(cv, "prior")$name), 
                 model = "cv.txt", 
                 n.iter = 11000, n.burnin = 1000, 
                 n.chains = 1, DIC = FALSE)
  
  # save mcmc results in appropiate directory
  param.mcmc <- as.mcmc(cv.bug$sims.matrix[, nodes(cv, "prior")$name])
  write.csv(param.mcmc, paste0("./AR", p, "/ew2007_param.csv"), row.names = FALSE)
  
  message(p)
}

##
## Calculate CV model density, q1(w_i), multivariate normal denstiy, q2(w_i), and their ratio, l(w_i), for both the bugs MCMC (w_1) and a simulated 
## mulitvariate normal sample (w_2)
##
library("tsbridge")

for (p in 0:8) {
  bugs.param <- read.csv(paste0("./AR", p, "/ew2007_param.csv"))
  
  # select random nodes
  cv <- ar.bugs(y = y, ar.order = p, beg = 9)
  theta <- subset(nodes(cv, part = "prior"), stoc == 1)$name
  bugs.param <- bugs.param[, theta, drop = FALSE]
  
  # adjust parameters to get onto continuous scale (-inf to +inf), then summary stats for MVN sample (w2)
  adju.param <- rescale(bug = cv, sims = bugs.param, to.real = TRUE)
  bugs.MU <- apply(adju.param, 2, mean)
  bugs.COV <- cov(adju.param)
  
  # take sample second sample (w2) and re-adjust to get back to same scale as bugs.param
  rmvn.param <- rmvnorm(10000, mean = bugs.MU, sigma = bugs.COV, method = "svd")
  rmvn.param <- rescale(bug = cv, sims = rmvn.param, to.real = FALSE)
  
  # calculate the mean (and variance of sv) at each time
  bugs.ymean <- y.fit(bug = cv, sims = bugs.param)
  rmvn.ymean <- y.fit(bug = cv, sims = rmvn.param)
  
  # redo bugs summary stats as previous were on adjusted parameters
  bugs.MU <- apply(bugs.param, 2, mean)
  bugs.COV <- cov(bugs.param)
  
  # calculate q1(w_i), q2(w_i) and l(w_i) for both the bugs mcmc (w_1) and mvn sample (w_2)
  w1 <- q1q2l(bug = cv, sims = bugs.param, ymean = bugs.ymean, MU = bugs.MU, COV = bugs.COV)
  w2 <- q1q2l(bug = cv, sims = rmvn.param, ymean = rmvn.ymean, MU = bugs.MU, COV = bugs.COV)
  
  # save information on w1 and w2 in appropriate model folder
  write.csv(w1, paste0("./AR", p, "/ew2007_w1.csv"), row.names = FALSE)
  write.csv(w2, paste0("./AR", p, "/ew2007_w2.csv"), row.names = FALSE)
  
  message(p)
}

##
## Calculate normalising constant and posterior model probability
##
nc <- data.frame(p = 0:8, r = NA, prob = NA)

# normalising constant
for (p in 0:8) {
  # read information on w1 and w2 from appropriate folder
  w1 <- read.csv(paste0("./AR", p, "/ew2007_w1.csv"))
  w2 <- read.csv(paste0("./AR", p, "/ew2007_w1.csv"))
  
  # calculate the normalising constant
  nc$r[p + 1] <- bridge(w1, w2, r0 = 500, tol = 1e-04)
  message(p)
}

# posterior model probability
nc.star <- max(nc$r)
nc$prob <- exp(nc$r - nc.star)/sum(exp(nc$r - nc.star))

# result
nc 

##
## Tidy up. Delete created folders and BUGS model
##
for (p in 0:8) {
  unlink(paste0("./AR", p), recursive = TRUE)
}
unlink("./cv.txt")