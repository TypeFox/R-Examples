############################################################################
#     MLwiN MCMC Manual
#
# 1   Introduction to MCMC Estimation and Bayesian Modelling . . . . . . . 1
#
#     Browne, W.J. (2009) MCMC Estimation in MLwiN, v2.13. Centre for
#     Multilevel Modelling, University of Bristol.
############################################################################
#     R script to replicate all analyses using R2MLwiN
#
#     Zhang, Z., Charlton, C., Parker, R, Leckie, G., and Browne, W.J.
#     Centre for Multilevel Modelling, 2012
#     http://www.bristol.ac.uk/cmm/software/R2MLwiN/
############################################################################

# 1.1 Bayesian modelling using Markov Chain Monte Carlo methods . . . . . .1

# 1.2 MCMC methods and Bayesian modelling . . . . . . . . . . . . . . . . .2

# 1.3 Default prior distributions . . . . . . . . . . . . . . . . . . . . .4

# 1.4 MCMC estimation . . . . . . . . . . . . . . . . . . . . . . . . . . .5

# 1.5 Gibbs sampling . . . . . . . . . . . . . . . . . . . . . . . . . . . 5

# 1.6 Metropolis Hastings sampling . . . . . . . . . . . . . . . . . . . . 8

# 1.7 Running macros to perform Gibbs sampling and Metropolis
#     Hastings sampling on the simple linear regression model . . . . . . 10

library(R2MLwiN)
# MLwiN folder
mlwin <- getOption("MLwiN_path")
while (!file.access(mlwin, mode = 1) == 0) {
  cat("Please specify the root MLwiN folder or the full path to the MLwiN executable:\n")
  mlwin <- scan(what = character(0), sep = "\n")
  mlwin <- gsub("\\", "/", mlwin, fixed = TRUE)
}
options(MLwiN_path = mlwin)


## Read tutorial data
data(tutorial, package = "R2MLwiN")

set.seed(1)

## Set variables
y <- tutorial$normexam
x <- tutorial$standlrt
N <- length(y)
xsq <- x^2
xy <- x * y
sumy <- sum(y)
sumx <- sum(x)
sumxsq <- sum(xsq)
sumxy <- sum(xy)

## Starting values for parameter estimates
beta0 <- 0
beta1 <- 0
sigma2e <- 1
epsilon <- 0.001
burnin <- 500
chain <- 5000
totaliterations <- burnin + chain
thinning <- 1
estimates <- matrix(, nrow = floor(chain/thinning), 3)
rownames(estimates) <- which((1:chain)%%thinning == 0)
colnames(estimates) <- c("beta0", "beta1", "sigma2e")

j <- 1
for (i in 1:totaliterations) {
  beta0 <- rnorm(1, (sumy - beta1 * sumx)/N, sqrt(sigma2e/N))
  beta1 <- rnorm(1, (sumxy - beta0 * sumx)/sumxsq, sqrt(sigma2e/sumxsq))
  
  e2i <- (y - (beta0 + beta1 * x))^2
  sume2i <- sum(e2i)
  
  sigma2e <- 1/rgamma(1, epsilon + N/2, epsilon + sume2i/2)
  
  if ((i%%thinning == 0) && (i > burnin)) {
    estimates[j, ] <- round(c(beta0, beta1, sigma2e), 3)
    j <- j + 1
  }
}

sumstat <- round(rbind(colMeans(estimates), apply(estimates, 2, sd)), 4)
rownames(sumstat) <- c("mean", "sd")
print(sumstat)

# 1.8 Dynamic traces for MCMC . . . . . . . . . . . . . . . . . . . . . . 12

par(mfrow = c(3, 1), mar = c(4, 4.5, 2, 2))
plot(1:nrow(estimates), estimates[, "beta0"], xlab = "iteration", ylab = expression(paste("Est. of ", beta[0])), type = "l")
plot(1:nrow(estimates), estimates[, "beta1"], xlab = "iteration", ylab = expression(paste("Est. of ", beta[1])), type = "l")
plot(1:nrow(estimates), estimates[, "sigma2e"], xlab = "iteration", ylab = expression(paste("Est. of ", sigma[e]^2)), 
  type = "l")

# 1.9 Macro to run a hybrid Metropolis and Gibbs sampling method
#     for a linear regression example . . . . . . . . . . . . . . . . . . 15
set.seed(1)

## Set variables
y <- tutorial$normexam
x <- tutorial$standlrt
N <- length(y)
xsq <- x^2
xy <- x * y
sumy <- sum(y)
sumx <- sum(x)
sumxsq <- sum(xsq)
sumxy <- sum(xy)

## Starting values for paramter estimates
beta0 <- 0
beta1 <- 0
sigma2e <- 1
epsilon <- 0.001
burnin <- 500
chain <- 5000
beta0sd <- 0.01
beta1sd <- 0.01
beta0accept <- 0
beta1accept <- 0
totaliterations <- burnin + chain
thinning <- 1
estimates <- matrix(, nrow = floor(chain/thinning), 3)
rownames(estimates) <- which((1:chain)%%thinning == 0)
colnames(estimates) <- c("beta0", "beta1", "sigma2e")

j <- 1
for (i in 1:totaliterations) {
  # Update beta0 Propose a new beta0
  beta0prop <- rnorm(1, beta0, beta0sd)
  beta0logpostdiff <- -1 * (2 * (beta0 - beta0prop) * (sumy - beta1 * sumx) + N * (beta0prop^2 - beta0^2))/(2 * 
    sigma2e)
  if (beta0logpostdiff > 0) {
    # Definitely accept as higher posterior
    beta0 <- beta0prop
    beta0accept <- beta0accept + 1
  } else {
    # Only sometimes accept
    if (runif(1) < exp(beta0logpostdiff)) {
      beta0 <- beta0prop
      beta0accept <- beta0accept + 1
    }
  }
  
  # Update beta1
  beta1prop <- rnorm(1, beta1, beta1sd)
  beta1logpostdiff <- -1 * (2 * (beta1 - beta1prop) * (sumxy - beta0 * sumx) + sumxsq * (beta1prop^2 - beta1^2))/(2 * 
    sigma2e)
  if (beta1logpostdiff > 0) {
    # Definitely accept as higher posterior
    beta1 <- beta1prop
    beta1accept <- beta1accept + 1
  } else {
    # Only sometimes accept
    if (runif(1) < exp(beta1logpostdiff)) {
      beta1 <- beta1prop
      beta1accept <- beta1accept + 1
    }
  }
  
  e2i <- (y - (beta0 + beta1 * x))^2
  sume2i <- sum(e2i)
  
  sigma2e <- 1/rgamma(1, epsilon + N/2, epsilon + sume2i/2)
  
  if ((i%%thinning == 0) && (i > burnin)) {
    estimates[j, ] <- round(c(beta0, beta1, sigma2e), 3)
    j <- j + 1
  }
}

sumstat <- round(rbind(colMeans(estimates), apply(estimates, 2, sd)), 4)
rownames(sumstat) <- c("mean", "sd")
print(sumstat)
cat(paste("The acceptance rate of beta0: ", round(beta0accept/totaliterations, 3), "\n"))
cat(paste("The acceptance rate of beta1: ", round(beta1accept/totaliterations, 3), "\n"))

# 1.10 MCMC estimation of multilevel models in MLwiN . . . . . . . . . . .18

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . . 19





############################################################################
