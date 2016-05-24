### R code from vignette source 'entropy_np.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: entropy_np.Rnw:70-72
###################################################
library(np)
options(prompt = "R> ", np.messages = FALSE, digits = 3)


###################################################
### code chunk number 2: entropy_np.Rnw:180-195
###################################################
set.seed(1234)
n <- 250
## Distributions are equal

sample.A <- data.frame(a=rnorm(n),b=factor(rbinom(n,2,.5)))
sample.B <- data.frame(a=rnorm(n),b=factor(rbinom(n,2,.5)))

npdeneqtest(sample.A,sample.B,boot.num=99)

## Distributions are unequal

sample.A <- data.frame(a=rnorm(n),b=factor(rbinom(n,2,.5)))
sample.B <- data.frame(a=rnorm(n,sd=10),b=factor(rbinom(n,2,.25)))

npdeneqtest(sample.A,sample.B,boot.num=99)


###################################################
### code chunk number 3: entropy_np.Rnw:227-255
###################################################
set.seed(1234)
n <- 1000
     
## Compute the statistic only, different distributions
     
x <- rchisq(n,df=10)
y <- rnorm(n,sd=10000,mean=-10000)
     
npunitest(x,y,bootstrap=FALSE)

## Data drawn from same continuous distribution

x <- rnorm(n)
y <- rnorm(n)
npunitest(x,y,boot.num=99)

## Data drawn from different continuous distributions having the
## same mean and variance

x <- rchisq(n,df=5)
y <- rnorm(n,mean=5,sd=sqrt(10))
npunitest(x,y,boot.num=99)

## Data drawn from different discrete distributions
     
x <- factor(rbinom(n,2,.5))
y <- factor(rbinom(n,2,.1))
npunitest(x,y,boot.num=99)


###################################################
### code chunk number 4: entropy_np.Rnw:313-326
###################################################
set.seed(1234)
     
n <- 100

## Asymmetric discrete probability distribution

x <- factor(rbinom(n,2,.8))
npsymtest(x,boot.num=99)

## Symmetric continuous distribution
     
y <- rnorm(n)
npsymtest(y,boot.num=99)


###################################################
### code chunk number 5: entropy_np.Rnw:366-375
###################################################
set.seed(123)
## Test/measure lack of fit between y and its fitted value from a
## regression model when x is relevant.
n <- 100
x <- rnorm(n)
y <- 1 + x + rnorm(n)
model <- lm(y~x)
y.fit <- fitted(model)
npdeptest(y,y.fit,boot.num=99,method="summation")


###################################################
### code chunk number 6: entropy_np.Rnw:409-424
###################################################
set.seed(123)
## A function to create a time series
ar.series <- function(phi,epsilon) {
  n <- length(epsilon)
  series <- numeric(n)
  series[1] <- epsilon[1]/(1-phi)
  for(i in 2:n) {
    series[i] <- phi*series[i-1] + epsilon[i]
  }
  return(series)
}
n <- 100
## Stationary persistent time-series
yt <- ar.series(0.95,rnorm(n))
npsdeptest(yt,lag.num=2,boot.num=99,method="summation")


###################################################
### code chunk number 7: entropy_np.Rnw:448-478
###################################################
Srho <- function(x,y,...) {
  ## First write a function to compute the integrand (this is fed to
  ## the `integrate' function). This function's first argument is the
  ## point at which the two densities are computed (the remaining
  ## arguments are the data vectors).
  integrand <- function(t,x,y) {
    ## First, nonparametrically estimate the density of the x data
    ## using a plug-in bandwidth and evaluate the density at the point
    ## `t'.
    f.x <- fitted(npudens(tdat=x,edat=t,bws=bw.SJ(x),...))
    ## Next, estimate the parametric density of the data y using the
    ## Gaussian distribution and evaluate the density at the point
    ## `t'.
    f.y <- dnorm(t,mean=mean(y),sd=sd(y))
    ## Compute and return the integrand evaluated at the point `t'.
    return(0.5*(sqrt(f.x)-sqrt(f.y))**2)
  }
  ## Feed the integrand function to integrate() and return the value
  ## of the integral.
  return(integrate(integrand,-Inf,Inf,x=x,y=y)$value)
}
set.seed(123)
n <- 1000
## Data drawn from the same distribution
x <- rnorm(n)
y <- rnorm(n)
Srho(x,y)
## Data drawn from different distributions
y <- rnorm(n,sd=100)
Srho(x,y)


