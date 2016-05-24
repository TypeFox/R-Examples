library(prinsimp)
library(MASS)

##  constructs  n random functions of the form
##     alpha0 + alpha1 * sin(2*pi*x/L) + alpha2 * sin(2*pi*x/M)
##            + sum_{i=1}^N  betai  sin( 2*pi*i*x/K)
##    with alpha0, alpha1, alpha2, betai independent normal random variables with mean 0 and
##            sds  a0,a1,a2,sd.sigma, e
##    evaluated at 1:x.max
##    add iid normal noise to each "data point"  with noise sd = e
##  output is a matrix of simulated data of dimension n by length(x.max)  
periodic.example <- function(L=72, M=5, K=12, a0, a1, a2, sd.sigma, e, n=100, x.max) {
    N <- length(sd.sigma)
    set.seed(123)
#
    beta1 <- rnorm(n*N, mean=rep(0,n*N), sd=rep(sd.sigma,n))
    beta <- t(matrix(beta1,nrow=N, ncol=n))
    
    ## alpha0 is coefficient for the constant
    alpha0 <- rnorm(n=n, mean=0, sd=a0)
    
    alpha1 <- rnorm(n=n, mean=0, sd=a1)
    
    alpha2 <- rnorm(n=n, mean=0, sd=a2)
    
    ## matrix of errors
    noise <- rnorm(n * x.max, mean=0, sd=e)
    error <- matrix(noise, nrow=n, ncol=x.max)
    
    ## construct the N repetitive parts
    x <- seq(1, x.max)
    x.argument <- (2 * pi * (1:N) %o% x) / K
    big.sin <- sin(x.argument)
    repetitive.parts <- beta %*% big.sin
    
    main.part <- alpha1 %o% sin(2*pi*x / L)
    
    other.part <- alpha2 %o% sin(2*pi*x / M)
    
    constant <- matrix(rep(alpha0, x.max), nrow=n, ncol=x.max)
    
    constant + main.part + other.part + repetitive.parts + error
}


## Simulates 100 curves of length 72 with the parameters in vignette
example <- periodic.example(a0=4, a1=4, a2=0, sd.sigma=.2, e=1, x.max=72)

## Plots the first 15 curves to see the periodic structure
matplot(t(example[1:15,]), type="l", xlab="x", ylab="y")

periodic.sim <- simpart(example, simpledim=70, "periodic", period=12)

## Analysis of the simulated periodic data with model space of
## dimension 2, nearly null space of dimension 70 with periodic
## simplicity measure and a period of 12.The two model basis vectors
## are plotted along with the variance-simplicity view and
## percent variance explained panel in a 2x2 grid with the layout argument
plot(periodic.sim, display=list(model=TRUE), layout=matrix(1:4, 2, 2))

## Demonstrate the use of na.action by deleting some elements
## of our data matrix
example.NA <- example
example.NA[2,3] <- NA

example.NA[53,30] <- NA

## Run simpart on the data matrix removing the missing data
periodic.omit <- simpart(example.NA, simpledim=70, "periodic", period=12, na.action=na.omit)

