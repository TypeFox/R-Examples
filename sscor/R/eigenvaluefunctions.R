## Calculates eigenvalues of SSCM;
# input:
# evShape: eigenvalues of shape matrix
# output: eigenvalues of SSCM

evShape2evSSCM <- function(evShape) {
  evShape <- evShape/sum(evShape)*length(evShape)
  delta <- fuint(evShape)*evShape/2
  return(delta)
}



## Calculates the eigenvalues of the shape matrix given the eigenvalues of the SSCM
# input: eigenvalues of the SSCM
# tol: absolute tolerance for the approximated eigenvalues
# itermax: maximal number of iterations of the approximation algorithm 
# output: list containing:
# eigenv: eigenvalues of the shape matrix
# iternumber: number of iterations of the approximation algorithm


evSSCM2evShape <- function(delta,tol=10^(-10),itermax=100) {
pm0 <- length(delta)
lambda <- numeric(pm0)
Index <- delta==0
p0 <- sum(Index)
delta <- delta[!Index]
p <- pm0-p0
if (p==2) {
lambda[1] <- 1
lambda[2] <- ((1-delta[1])/(delta[1]))^2
return(lambda)
}

pdata <- p
pdata <- min(pdata,200)
jacobiquad <- get(load(system.file("extdata", "jacobiquad", package = "sscor")))	
kno <- jacobiquad[[1]][pdata,]
wei <- jacobiquad[[2]][pdata,]

lambdaalt <- lambdaneu <- delta
lambdaalt[1] <- lambdaalt[1]+1

for(i in 1:itermax) {
lambdaalt <- lambdaneu
integr <- fuint(lambdaalt,wei=wei,kno=kno)
lambdaneu <- 2*delta/integr
lambdaneu <- lambdaneu/sum(lambdaneu)*p
deltaapprox <- lambdaalt*integr/2
if(sum(abs(deltaapprox-delta))<tol) break
if(i==itermax) warning("Maximal number of iterations reached. Approximated shape matrix might be imprecise.")
}
lambda[!Index] <- lambdaneu/sum(lambdaneu)
return(lambda)
}

## Evaluates the integral which gives the eigenvalues of the SSCM;
# input:
# evShape: eigenvalues of shape matrix
# wei: weights of Gauss-Jacobi quadrature
# kno: knots of Gauss-Jacobi quadrature
# output: defined integral


fuint <- function(lambda,wei=NULL,kno=NULL) {
pdata <- length(lambda)
pdata <- min(pdata,200)
p <- pdata/2+1
if(is.null(wei)|is.null(kno)) {
jacobiquad <- get(load(system.file("extdata", "jacobiquad", package = "sscor")))
kno <- jacobiquad[[1]][pdata,]
wei <- jacobiquad[[2]][pdata,]
}
wei <- wei*(1-kno)^(-p)
kno <- (1+kno)/(1-kno)
foo <- function(x){ 1/(sqrt(apply(1+outer(lambda,x),2,prod)))}
f1 <- foo(kno)
fd <- 1/(1+outer(lambda,kno))
w <- f1*t(fd)
bi <- 2*apply(wei*w,2,sum)
return(bi)
}

