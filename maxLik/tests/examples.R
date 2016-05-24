library( maxLik )
options(digits=4)

printRounded <- function( x ) {
   for( i in names( x ) ) {
      cat ( "$", i, "\n", sep = "" )
      if( is.numeric( x[[i]] ) ) {
         print( round( x[[i]], 3 ) )
      } else {
         print( x[[i]] )
      }
      cat( "\n" )
   }
   cat( "attr(,\"class\")\n" )
   print( class( x ) )
}

# round gradients to increase reproducibility of the accuracy
roundGradients <- function( object ) {
   object$gradient <- round( object$gradient, 3 )
   return( object )   
}

### activePar
# a simple two-dimensional exponential hat
f <- function(a) exp(-(a[1]-2)^2 - (a[2]-4)^2)
#
# maximize wrt. both parameters 
free <- maxNR(f, start=1:2)
printRounded( free )
free <- roundGradients( free )
summary(free)  # results should be close to (2,4)
activePar(free)
# allow only the second parameter to vary
cons <- maxNR(f, start=1:2, activePar=c(FALSE,TRUE))
printRounded( cons )
cons <- roundGradients( cons )
summary(cons) # result should be around (1,4)
activePar(cons)
# specify fixed par in different ways
cons2 <- maxNR(f, start=1:2, fixed=1)
all.equal( cons[-3], cons2[-3] )
cons3 <- maxNR(f, start=1:2, fixed=c(TRUE,FALSE))
all.equal( cons[-3], cons3[-3] )
cons4 <- maxNR(f, start=c(a=1, b=2), fixed="a")
cons4 <- roundGradients( cons4 )
print(summary(cons4))
all.equal( cons, cons4 )

### compareDerivatives
set.seed( 2 )
## A simple example with sin(x)' = cos(x)
f <- sin
compareDerivatives(f, cos, t0=1)
##
## Example of log-likelihood of normal density.  Two-parameter
## function.
x <- rnorm(100, 1, 2) # generate rnorm x
l <- function(b) sum(log(dnorm((x-b[1])/b[2])/b[2]))
              # b[1] - mu, b[2] - sigma
gradl <- function(b) {
   c(sum(x - b[1])/b[2]^2,
   sum((x - b[1])^2/b[2]^3 - 1/b[2]))
}
compareDerivatives(l, gradl, t0=c(1,2))


### hessian
set.seed( 3 )
# log-likelihood for normal density
# a[1] - mean
# a[2] - standard deviation
ll <- function(a) sum(-log(a[2]) - (x - a[1])^2/(2*a[2]^2))
x <- rnorm(1000) # sample from standard normal
ml <- maxLik(ll, start=c(1,1))
# ignore eventual warnings "NaNs produced in: log(x)"
printRounded( ml )
print( ml )
summary(ml) # result should be close to c(0,1)
hessian(ml) # How the Hessian looks like
sqrt(-solve(hessian(ml))) # Note: standard deviations are on the diagonal
print(stdEr(ml))
                           # test vector of stdEr-s
#
# Now run the same example while fixing a[2] = 1
mlf <- maxLik(ll, start=c(1,1), activePar=c(TRUE, FALSE))
printRounded( mlf )
print( mlf )
summary(mlf) # first parameter close to 0, the second exactly 1.0
hessian(mlf)
# now invert only the free parameter part of the Hessian
sqrt(-solve(hessian(mlf)[activePar(mlf), activePar(mlf)]))
# gives the standard deviation for the mean
print(stdEr(mlf))
                           # test standard errors with fixed par


### maxBFGS
set.seed( 5 )
# Maximum Likelihood estimation of the parameter of Poissonian distribution
n <- rpois(100, 3)
loglik <- function(l) n*log(l) - l - lfactorial(n)
# we use numeric gradient
a <- maxBFGS(loglik, start=1)
print( a[-3] )
class( a )
a <- roundGradients( a )
summary( a )
# you would probably prefer mean(n) instead of that ;-)
# Note also that maxLik is better suited for Maximum Likelihood


### logLik.maxLik
set.seed( 4 )
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
gradlik <- function(theta) 1/theta - t
hesslik <- function(theta) -100/theta^2
## Estimate with analytic gradient and hessian
a <- maxLik(loglik, gradlik, hesslik, start=1)
printRounded( a )
print( a )
## print log likelihood value
logLik( a )
## compare with log likelihood value of summary object
all.equal( logLik( a ), logLik( summary( a ) ) )


### maxBHHH
set.seed( 6 )
## ML estimation of exponential duration model:
t <- rexp(100, 2)
## Estimate with numeric gradient and hessian
a <- maxBHHH(loglik, start=1, print.level=2)
print( a )
a <- roundGradients( a )
summary(a)
## Estimate with analytic gradient
a <- maxBHHH(loglik, gradlik, start=1)
print( a )
a <- roundGradients( a )
summary(a)


### maxLik
set.seed( 7 )
## ML estimation of exponential duration model:
t <- rexp(100, 2)
## Estimate with numeric gradient and hessian
a <- maxLik(loglik, start=1, print.level=2)
printRounded( a )
print( a )
summary(a)
## Estimate with analytic gradient and hessian
a <- maxLik(loglik, gradlik, hesslik, start=1)
printRounded( a )
print( a )
summary(a)


### maxNR
set.seed( 8 )
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglikSum <- function(theta) sum(log(theta) - theta*t)
## Note the log-likelihood and gradient are summed over observations
gradlikSum <- function(theta) sum(1/theta - t)
## Estimate with numeric gradient and Hessian
a <- maxNR(loglikSum, start=1, print.level=2)
a <- roundGradients( a )
print( a )
summary(a)
## You would probably prefer 1/mean(t) instead ;-)
## Estimate with analytic gradient and Hessian
a <- maxNR(loglikSum, gradlikSum, hesslik, start=1)
a <- roundGradients( a )
print( a )
summary(a)


### maximType
## maximise two-dimensional exponential hat.  Maximum is at c(2,1):
f <- function(a) exp(-(a[1] - 2)^2 - (a[2] - 1)^2)
m <- maxNR(f, start=c(0,0))
m <- roundGradients( m )
print( m )
summary(m)
maximType(m)
## Now use BFGS maximisation.
m <- maxBFGS(f, start=c(0,0))
m <- roundGradients( m )
print( m )
summary(m)
maximType(m)

### Test maxNR with 0 iterations.  Should perform no iterations
### Request by Yves Croissant
f <- function(a) exp(-(a[1] - 2)^2 - (a[2] - 1)^2)
m0 <- maxNR(f, start=c(1.1, 2.1), iterlim=0)
m0 <- roundGradients( m0 )
summary(m0)

### nObs
set.seed( 10 )
# Construct a simple OLS regression:
x1 <- runif(100)
x2 <- runif(100)
y <- 3 + 4*x1 + 5*x2 + rnorm(100)
m <- lm(y~x1+x2)  # estimate it
nObs(m)


### nParam
set.seed( 11 )
# Construct a simple OLS regression:
x1 <- runif(100)
x2 <- runif(100)
y <- 3 + 4*x1 + 5*x2 + rnorm(100)
m <- lm(y~x1+x2)  # estimate it
summary(m)
nParam(m) # you get 3


### numericGradient
# A simple example with Gaussian bell
f0 <- function(t0) exp(-t0[1]^2 - t0[2]^2)
numericGradient(f0, c(1,2))
numericHessian(f0, t0=c(1,2))
# An example with the analytic gradient
gradf0 <- function(t0) -2*t0*f0(t0)
numericHessian(f0, gradf0, t0=c(1,2))
# The results should be similar as in the previous case
# The central numeric derivatives have usually quite a high precision
compareDerivatives(f0, gradf0, t0=1:2)
# The differenc is around 1e-10


### returnCode
## maximise the exponential bell
f1 <- function(x) exp(-x^2)
a <- maxNR(f1, start=2)
a <- roundGradients( a )
print( a )
returnCode(a) # should be success (1 or 2)
## Now try to maximise log() function
a <- maxNR(log, start=2)
returnCode(a) # should give a failure (4)


### returnMessage
## maximise the exponential bell
f1 <- function(x) exp(-x^2)
a <- maxNR(f1, start=2)
a <- roundGradients( a )
print( a )
returnMessage(a) # should be success (1 or 2)
## Now try to maximise log() function
f2 <- function(x) log(x)
a <- maxNR(log, start=2)
returnMessage(a) # should give 'Iteration limit exceeded.'


### summary.maxLik
set.seed( 15 )
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
gradlik <- function(theta) 1/theta - t
hesslik <- function(theta) -100/theta^2
## Estimate with numeric gradient and hessian
a <- maxLik(loglik, start=1, print.level=2)
printRounded( a )
print( a )
summary(a)
## Estimate with analytic gradient and hessian
a <- maxLik(loglik, gradlik, hesslik, start=1)
printRounded( a )
print( a )
summary(a)


### summary.maxim and for "gradient"/"hessian" attributes
### Test for infinity
## maximize a 2D quadratic function:
f <- function(b) {
  x <- b[1]; y <- b[2];
    val <- (x - 2)^2 + (y - 3)^2
    attr(val, "gradient") <- c(2*x - 4, 2*y - 6)
    attr(val, "hessian") <- matrix(c(2, 0, 0, 2), 2, 2)
    val
}
## Use c(0,0) as initial value.  
result1 <- maxNR( f, start = c(0,0) )
print( result1 )
summary( result1 )
## Now use c(1000000, -777777) as initial value and ask for hessian
result2 <- maxNR( f, start = c( 1000000, -777777))
print( result2 )
summary( result2 )


### Test for "gradient"/"hessian" attributes.  A case which converges.
hub <- function(x) {
   v <- exp(-sum(x*x))
   val <- v
   attr(val, "gradient") <- -2*x*v
   attr(val, "hessian") <- 4*(x %*% t(x))*v - diag(2*c(v, v))
   val
}
a <- maxNR(hub, start=c(2,1))
a <- roundGradients( a )
summary( a )
## Now test "gradient" attribute for BHHH/3-parameter probit
N <- 1000
loglikProbit <- function( beta) {
   xb <- x %*% beta
   loglik <- ifelse(y == 0,
                    pnorm( xb, log=TRUE, lower.tail=FALSE),
                    pnorm( xb, log.p=TRUE))
   grad <- ifelse(y == 0,
                  -dnorm(xb)/pnorm(xb, lower.tail=FALSE),
                  dnorm(xb)/pnorm(xb))
   grad <- grad*x
   attr(loglik, "gradient") <- grad
   loglik
}
x <- runif(N)
x <- cbind(x, x - runif(N), x - runif(N))
y <- x[,1] + 2*x[,2] - x[,3] + rnorm(N) > 0
summary(maxLik(loglikProbit, start=c(0,0,0), method="bhhh"))



### vcov.maxLik
set.seed( 17 )
## ML estimation of exponential duration model:
t <- rexp(100, 2)
## Estimate with numeric gradient and hessian
a <- maxLik(loglik, start=1, print.level=2)
printRounded( a )
print( a )
round( vcov( a ), 3 )
## Estimate with analytic gradient and hessian
a <- maxLik(loglik, gradlik, hesslik, start=1)
printRounded( a )
print( a )
round( vcov( a ), 3 )
print(stdEr(a))
                           # test single stdEr
