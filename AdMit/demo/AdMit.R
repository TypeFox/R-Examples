wait <- function()
{
  t <- readline("\nPlease 'q' to quit the demo or any other key to continue...\n")
  if (t == "q") stop ("end of the demo")
}

#######################################################################
## START OF THE DEMO
#######################################################################

##_____________________________________________________________________
## GELMAN AND MENG (1991) EXAMPLE

## Initialization
options(digits = 4, max.print = 40, prompt = "R> ")

## Gelman and Meng (1991) kernel function
GelmanMeng <- function(x, A = 1, B = 0, C1 = 3, C2 = 3, log = TRUE)
{
  if (is.vector(x))
    x <- matrix(x, nrow = 1)
  r <- -.5 * (A * x[,1]^2 * x[,2]^2 + x[,1]^2 + x[,2]^2
            - 2 * B * x[,1] * x[,2] - 2 * C1 * x[,1] - 2 * C2 * x[,2])
  if (!log)
    r <- exp(r)
  as.vector(r)
}
wait()

## Contour plot of the Gelman and Meng (1991) kernel function.
PlotGelmanMeng <- function(x1, x2)
{
  GelmanMeng(cbind(x1, x2), log = FALSE)
}
x1 <- x2 <- seq(from = -1.0, to = 6.0, by = 0.1)
z <- outer(x1, x2, FUN = PlotGelmanMeng)
contour(x1, x2, z, nlevel = 20, las = 1, lwd = 2,
        col = rainbow(20),
        xlab = expression(X[1]), ylab = expression(X[2]))
abline(a = 0, b = 1, lty = "dotted")
wait()

## Use AdMit to the Gelman and Meng (1991) kernel function.
set.seed(1234)
outAdMit <- AdMit(GelmanMeng, mu0 = c(0.0, 0.1))
print(outAdMit)
wait()

## Contour plot of the mixture approximation obtained with AdMit.
PlotMit <- function(x1, x2, mit)
{
  dMit(cbind(x1, x2), mit = mit, log = FALSE)
}
z <- outer(x1, x2, FUN = PlotMit, mit = outAdMit$mit)
contour(x1, x2, z, nlevel = 20, las = 1, lwd = 2, 
        col = rainbow(20),
        xlab = expression(X[1]), ylab = expression(X[2]))
abline(a = 0, b = 1, lty = "dotted")
wait()

## Contour plot of the components of the mixture approximation
## obtained with AdMit.
par(mfrow = c(2,2))
for (h in 1:4)
{
  mith <- list(p = 1,
               mu = outAdMit$mit$mu[h,,drop = FALSE],
               Sigma = outAdMit$mit$Sigma[h,,drop = FALSE],
               df = outAdMit$mit$df)
  z <- outer(x1, x2, FUN = PlotMit, mit = mith)
  contour(x1, x2, z, las = 1, nlevel = 20, lwd = 2,
          col = rainbow(20),
          xlab = expression(X[1]), ylab = expression(X[2]))
    abline(a = 0, b = 1, lty = "dotted")
    title(main = paste("component nr.", h))
}
wait()

## Use importance sampling with the mixture approximation
## as the importance density.
outAdMitIS <- AdMitIS(KERNEL = GelmanMeng, mit = outAdMit$mit)
print(outAdMitIS)
wait()

## Use an alternative 'G' function in importance sampling
## for computing the covariance matrix.
G.cov <- function(theta, mu)
{
  G.cov_sub <- function(x)
    (x-mu) %*% t(x-mu)

  theta <- as.matrix(theta)
  tmp <- apply(theta, 1, G.cov_sub)
  if (length(mu) > 1)
    t(tmp)
  else
    as.matrix(tmp)
}
outAdMitIS <- AdMitIS(KERNEL = GelmanMeng, G = G.cov, mit = outAdMit$mit, mu = c(1.459, 1.459))
print(outAdMitIS)
V <- matrix(outAdMitIS$ghat, 2, 2)
print(V)
cov2cor(V)
wait()

## Use independence Metropolis-Hastings algorithm with
## the mixture approximation as the candidate density.
outAdMitMH <- AdMitMH(KERNEL = GelmanMeng, mit = outAdMit$mit)
print(outAdMitMH)
wait()

## Use some functions of the package 'coda' to obtain summaries from the MCMC output.
library("coda")
draws <- as.mcmc(outAdMitMH$draws[1001:1e5,])
colnames(draws) <- c("X1", "X2")
summary(draws)$stat
summary(draws)$stat[,3]^2 / summary(draws)$stat[,4]^2
plot(draws)
wait()

##_____________________________________________________________________
## SIMPLE ECONOMETRIC EXAMPLE

## Simple econometric model: y_t ~ i.i.d. N(mu,sigma^2)
## Jeffreys prior p(theta) prop 1/sigma for sigma > 0
KERNEL <- function(theta, y, log = TRUE)
{
  if (is.vector(theta))
    theta <- matrix(theta, nrow = 1)

  KERNEL_sub <- function(thetai)
  {
    if (thetai[2] > 0)
    {
      r <- - log(thetai[2]) + sum(dnorm(y, thetai[1], thetai[2], TRUE))
    }
    else
    {
      r <- -Inf
    }
    as.numeric(r)
  }
  r <- apply(theta, 1, KERNEL_sub)
  if (!log)
    r <- exp(r)
  as.numeric(r)
}

## Generate 20 draws for mu = 1 and sigma = 0.5
set.seed(1234)
y <- rnorm(20, 2.0, 0.5)
par(mfrow = c(1,1))
plot(y)
wait()

## Run AdMit (with default values); pass the vector y
## of observations using the ... argument of AdMit and
## print steps of the fitting process
outAdMit <- AdMit(KERNEL, mu0 = c(1.0, 1.0), y = y, control = list(trace = TRUE))
print(outAdMit)
wait()

## Use independence Metropolis-Hastings algorithm with
## the mixture approximation as the candidate density; pass the
## vector y of observations using the ... argument of AdMitMH
outAdMitMH <- AdMitMH(KERNEL = KERNEL, mit = outAdMit$mit, y = y)
print(outAdMitMH)
wait()

## Use some functions of the package 'coda' to obtain summaries from the MCMC output.
draws <- as.mcmc(outAdMitMH$draws[1001:1e5,])
colnames(draws) <- c("mu", "sigma")
summary(draws)$stat
summary(draws)$stat[,3]^2 / summary(draws)$stat[,4]^2
plot(draws)
wait()

##_____________________________________________________________________
## BAYESIAN ESTIMATION OF A MIXTURE OF ARCH(1) MODEL

## Define the prior density
## The function outputs a Nx2 matrix. The first column indicates whether the
## prior constraint is satisfied, the second returns the value of the prior
PRIOR <- function(omega1, omega2, alpha, p, log = TRUE)
{
  c1 <- (omega1 > 0.0 & omega2 > 0.0 & alpha >= 0.0)   ## positivity constraint
  c2 <- (alpha < 1.0)                                  ## stationarity constraint
  c3 <- (p > 0.0 & p < 1.0)                            ## U(0,1) prior on p
  c4 <- (omega1 < omega2)                              ## identification constraint
  r1 <- c1 & c2 & c3 & c4
  r2 <- rep.int(-Inf, length(omega1))
  tmp <- dnorm(omega1[r1==TRUE], 0.0, 2.0, log = TRUE)       ## prior on omega1
  tmp <- tmp + dnorm(omega2[r1==TRUE], 0.0, 2.0, log = TRUE) ## prior on omega2
  r2[r1==TRUE] <- tmp + dnorm(alpha[r1==TRUE], 0.2, 0.5, log = TRUE) ## prior on alpha
  if (!log)
    r2 <- exp(r2)
  cbind(r1, r2)
}
wait()

## Define the kernel function
## The function takes a Nx4 matrix of parameters (theta), a vector of log-returns (y)
## It outputs the kernel value for the N parameters
KERNEL <- function(theta, y, log = TRUE)
{
  if (is.vector(theta))
    theta <- matrix(theta, nrow = 1)
  N <- nrow(theta)

  ## compute the prior for the parameters
  prior <- PRIOR(theta[,1], theta[,2], theta[,3], theta[,4])

  ## the kernel function is implemented in C in order to speed up the estimation
  d <- .C(name = "fnKernelMixtureArch_C",
          theta = as.double( as.vector(t(theta)) ),
          N = as.integer(N),
          y = as.double(y),
          n = as.integer(length(y)),
          prior = as.double( as.vector(t(prior)) ),
          d = vector("double", N),
          PACKAGE = "AdMit",
          NAOK = TRUE,
          DUP = FALSE)$d
    
    if (!log)
      d <- exp(d)
    as.vector(d)
  }
wait()

## Load the data set
library("fEcofin")
data("dem2gbp")
y <- dem2gbp[1:250,1]
par(mfrow = c(1,1))
plot(y, type = "l", las = 1, ylab = "log-returns", xlab = "time index")
wait()

## Maximize to find the mode of the kernel function
NLL <- function(..., log = TRUE) -KERNEL(...)
start <- c(0.1, 0.5, 0.1, 0.5)
outML <- optim(par = start, fn = NLL, y = y, method = "Nelder-Mead",
               control = list(trace = 1, maxit = 5000))
## print the mode
round(outML$par, 4)

## __Adaptive mixture approach___
set.seed(1234)
system.time( outAdMit <- AdMit(KERNEL, mu0 = outML$par, y = y, control = list(IS = TRUE, trace = TRUE)) )
print(outAdMit)

## ___Naive (unimodal Student-t approach) approach___
set.seed(1234)
system.time( outAdMit <- AdMit(KERNEL, mu0 = outML$par, y = y, control = list(Hmax = 1)) )
print(outAdMit)

## Then, use the output of 'AdMit' to perform importance sampling or independence MH sampling

## __Importance sampling approach (for estimating the posterior mean)__
set.seed(1234)
outAdMitIS <- AdMitIS(N = 50000, KERNEL = KERNEL, mit = outAdMit$mit, y = y)
print(outAdMitIS)

## __Importance sampling approach (for estimating the posterior covariance matrix)__
set.seed(1234)
## !!! compile the G.cov function above (section 3) !!!
outAdMitIS <- AdMitIS(N = 50000, KERNEL = KERNEL, G = G.cov, mit = outAdMit$mit, y = y, mu = outAdMitIS$ghat)
print(outAdMitIS)
## posterior standard deviations
sqrt(diag(matrix(outAdMitIS$ghat, 4, 4)))

## __Independence chain Metropolis-Hasting algorithm__
set.seed(1234)
outAdMitMH <- AdMitMH(N = 51000, KERNEL = KERNEL, mit = outAdMit$mit, y = y)
print(outAdMitMH$accept)
draws <- outAdMitMH$draws[1001:nrow(outAdMitMH$draws),]
colnames(draws) <- c("omega1", "omega2", "alpha", "p")

## ACF plots of the MCMC output
par(mfrow = c(2,2))
par(cex.axis = 1.2, cex.lab = 1.2)
acf(draws[,"omega1"], lag.max = 30, las = 1, main = expression(omega[1]))
acf(draws[,"omega2"], lag.max = 30, las = 1, main = expression(omega[2]))
acf(draws[,"alpha"], lag.max = 30, las = 1, main = expression(alpha))
acf(draws[,"p"], lag.max = 30, las = 1, main=expression(p))
## ACF up to lag 10
apply(draws, 2, acf, plot = FALSE, lag.max = 10)

## use summary from package coda
draws <- as.mcmc(draws)
summary(draws)$stat
## RNE
summary(draws)$stat[,3]^2 / summary(draws)$stat[,4]^2

#######################################################################
## END OF THE DEMO
## Additional code in ./AdMitJSS.R and ./AdMitRnews.R files
#######################################################################