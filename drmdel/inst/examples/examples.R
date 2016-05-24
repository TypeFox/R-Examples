# Load the drmdel package
library(drmdel)

# ##############################
# Data generation
set.seed(25)
n_samples <- c(100, 200, 180, 150, 175)  # sample sizes
x0 <- rgamma(n_samples[1], shape=5, rate=1.8)
x1 <- rgamma(n_samples[2], shape=12, rate=1.2)
x2 <- rgamma(n_samples[3], shape=12, rate=1.2)
x3 <- rgamma(n_samples[4], shape=18, rate=5)
x4 <- rgamma(n_samples[5], shape=25, rate=2.6)
x <- c(x0, x1, x2, x3, x4)
# ##############################

# ##############################
# Fit a DRM with the basis function q(x) = (x, log(abs(x))),
# which is the basis function for gamma family.

# There are 11 built-in basis function in drmdel(). And q(x)
# = (x, log(abs(x))) is the 6th basis function, so we can
# fit the model by specifying basis_func=6 in drmdel() as
# follows:
drmfit <- drmdel(x=x, n_samples=n_samples, basis_func=6)

# A brief summary of the DRM fit
summaryDRM(drmfit)

# List the names of the outputs
names(drmfit)

# Estimate the covariance matrix of the maximum empirical
# likelihood estimator (MELE)
meleCov(drmfit)

# Another way of specifying basis function for drmdel() is
# to pass a user-specified R function to the basis_func
# argument of the drmdel() function.
basis_gamma <- function(x) return(c(x, log(abs(x))))
drmfit1 <- drmdel(x=x, n_samples=n_samples, basis_func=basis_gamma)

# One can see the summary of this DRM fit is exactly the
# same as that of the previous fit with basis_func=6
summaryDRM(drmfit1)

# NOTE: If the basis function one wants to use is included
# in the built-in function list, one should use the built-in
# functions by passing an integer between 1 to 11 to the
# drmdel() function, because the computation will be faster
# with a built-in function.
# ##############################

# ##############################
# Testing the equality of all the distribution functions

# Suppose we want to test:
#     H_0: F_0=F_1=...=F_4
#         against
#     H_1: One distribution functin is not the same as the others

# The results of the dual empirical likelihood ratio (DELR)
# test for this hypothesis is already given in the output of
summaryDRM(drmfit)

# The output for DELR test looks like:
# Default hypothesis testing problem:
#     H_0: All distribution functions, F_k's, are equal.
#     H_1: One of the distribution functions is different from the others.
# Dual empirical likelihood ratio statistc (DELR): 1035.261 
# Degree of freedom: 8 
# p-value: 0

# NOTE: The drmdel() function can be used to test more
# complicated hypothesis about DRM parameter beta. We will
# illustrate that in an upcoming article about the DRM
# package.
# ##############################

# ##############################
# Quantile estimation

# We now have five samples labeled from 0 -- 4. Note that
# the base-line distribution is always labeled as 0.

# Denote the p^th quantile of the k^th, k=0, 1, ..., 4,
# distribution as q_{k,p}. To estimate q_{0,0.25},
# q_{0,0.6}, q_{1,0.1}, q_{2,0.1}, we use the following
# function
qe <- quantileDRM(k=c(0, 0, 1, 2), p=c(0.25, 0.6, 0.1, 0.1),
                      drmfit=drmfit)

# Print the output
qe
# There are there outputs:
#   est -- quantiles estimates
#
#   cov -- estimated covariance matrix of the quantile
#     estimators
# ##############################

# ##############################
# Quantile comparison

# ----------Example 1----------
# Suppose we want to test 
#     H_0: q_{1,0.05} = q_{2,0.05}
#           against
#     H_0: q_{1,0.05} != q_{2,0.05}
#
# It is known that
#     sqrt(n)*(q_{1,0.05} - qHat_{1,0.05}, qHat_{2,0.05} - q_{2,0.05}),
# where n is the total sample size, has a bivariate normal
# distribution. So under the null,
#     sqrt(n)*(qHat_{1,0.05} - qHat_{2,0.05})/sqrt({(1, -1)^T sigmaHat(1, -1)}),
# where sigmaHat is the estimated asymptotic covariance
# matrix of these two quantile estimators, has a asymptotic
# standard normal distribution. This gives us a Wald-type
# test for comparing two quantiles.

# Calculate quantile estimator and the corresponding
# covariance matrix
qe <- quantileDRM(k=c(1, 2), p=c(0.05, 0.05), drmfit=drmfit)
qe

n_total <- sum(n_samples)
# Define test-statistic
( wald_stat <- sqrt(n_total) * (qe$est[1] - qe$est[2])
    / sqrt( rbind(c(-1,1))%*% qe$cov %*% cbind(c(-1,1)) ) )

# p-value
(p_val <- 1 - pnorm(abs(wald_stat)))
# P-value is 0.269. Don't reject the null.
# --------------------

# ----------Example 2----------
# Suppose now we want to test 
#     H_0: q_{0,0.05} = q_{2,0.05}
#           against
#     H_0: q_{0,0.05} != q_{2,0.05}

# Calculate quantile estimator and the corresponding
# covariance matrix
qe1 <- quantileDRM(k=c(0, 2), p=c(0.05, 0.05), drmfit=drmfit)
qe1

# Define test-statistic
( wald_stat1 <- sqrt(n_total) * (qe1$est[1] - qe1$est[2])
    / sqrt( rbind(c(-1,1))%*% qe1$cov %*% cbind(c(-1,1)) ) )

# p-value
(p_val1 <- 1 - pnorm(abs(wald_stat1)))
# P-value is 0. Reject the null.
# ---------------
# ##############################

# ##############################
# Density estimation

# To estimate the density of population 3 under
# the DRM, we use
dens_pop3 <- densityDRM(k=3, drmfit=drmfit)
# The output is a an object with class '"density"'. See
# ?density for help.

# Plot the estimated density
plot(dens_pop3, main=bquote(F[3]), ylim=range(c(0, 0.5)))
# Don't forget to specify 'main' argument to plot(),
# otherwise you will see a super long plot title covering
# the whole plot!

# Add the empirical kernel density estimation curve of F_3
# based on x3 on the above density plot
lines(density(x3), col="green4", lty="28F8")

# Add the true density curve of F_3 on the above density
# plot
lines(seq(min(dens_pop3$y), max(dens_pop3$x), 0.01),
      dgamma(seq(min(dens_pop3$y), max(dens_pop3$x), 0.01), 18, 5),
      type="l", col="red", lty="dotted")

legend(9, 0.5,
       legend=c("DRM density estimator", "Empirical kernel density estimator",
                "True density"),
       col=c("black", "green4", "red"), lty=c("solid", "28F8", "dotted"))
# ##############################

# ##############################
# Cumulative distribution function (CDF) estimation

# To estimate the CDF of population 1 at c(3, 7.5, 11) and
# that of population 3 at c(2, 6), we use
cdf_est1 <- cdfDRM(k=c(1, 3), x=list(c(3, 7.5, 11), c(2, 6)),
                   drmfit=drmfit)
# list the output
names(cdf_est1)
cdf_est1$F1
cdf_est1$F3

# To estimate the CDF of population 2 and 4 at the observed
# data points under the DRM, we use
cdf_est <- cdfDRM(k=c(2, 4), drmfit=drmfit)
# list the output
names(cdf_est)
cdf_est$F2
cdf_est$F4
# ##############################
