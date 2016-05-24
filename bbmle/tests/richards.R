## implement richards-incidence (="revised superlogistic")
##  with analytic gradients

## from Junling's code:
model_richardson <- function(times, theta, N)
{	
	x0 = theta[1]
	lambda = theta[2]
	K = theta[3] * N
	alpha = theta[4]
	return(K/(1+((K/x0)^alpha-1)*exp(-lambda*alpha*times))^(1/alpha))
}

## equivalent model, in terms of sigma and as a symbolic expression
Rcum <- expression((sigma*N)/(1+(((sigma*N)/x0)^alpha-1)*exp(-lambda*alpha*times))^(1/alpha))

pnames <- c("x0","lambda","sigma","alpha")

## function to compute gradient (and value), derived by R
Rderiv <- deriv(Rcum,pnames, function.arg=c(pnames,"N","times"))

## equivalent (using Rcum): return incidence (incid=TRUE) or cumulative incidence (incid=FALSE)
calc_mean <- function(p,times,N,incid=TRUE) {
  ## this is more 'magic' than I would like it to be ...
  ##  have to create an environment and populate it with the contents of p (and N and times),
  ##  then evaluate the expression in this environment
  pp <- c(as.list(p),list(times=times,N=N))
  ##  e0 <- new.env()
  ## mapply(assign,names(pp),pp,MoreArgs=list(envir=e0))
  cumvals <- eval(Rcum,envir=pp)
  if (incid) diff(cumvals) else cumvals
}

## Poisson likelihood function
likfun <- function(p,dat,times,N,incid=TRUE) {
  -sum(dpois(dat,calc_mean(p,times,N,incid=incid),log=TRUE))
}

## deriv of P(x,lambda) = -sum(dpois(x,lambda,log=TRUE)) wrt lambda == sum(1-lambda/x) = N - lambda/(sum(x))
## deriv of P(x,lambda) wrt p = dP/d(lambda) * d(lambda)/dp

## compute gradient vector
gradlikfun <- function(p,dat,times,N,incid=TRUE) {
  gcall <- do.call(Rderiv,c(as.list(p),list(times=times,N=N))) ## values + gradient matrix
  lambda <- gcall
  attr(lambda,"gradient") <- NULL
  if (incid) lambda <- diff(lambda)
  gmat <- attr(gcall,"gradient") ## extract gradient
  if (incid) gmat <- apply(gmat,2,diff)  ## differences
  totderiv <- sweep(gmat,MARGIN=1,(1-dat/lambda),"*") ## apply chain rule (multiply columns of gmat by dP/dlambda)
  colSums(totderiv)  ## deriv of summed likelihood = sum of derivs of likelihod
}

N <- 1000
p0 <- c(x0=0.1,lambda=1,sigma=0.5,alpha=0.5)
t0 <- 1:10

## deterministic versions of data  (cumulative and incidence)
dcdat <- model_richardson(t0,p0,N)
ddat <- diff(dcdat)

plot(t0,dcdat)
plot(t0[-1],ddat)

set.seed(1001)
ddat <- rpois(length(ddat),ddat)

likfun(p0,ddat,t0,N)
gradlikfun(p0,ddat,t0,N)

library(numDeriv)
grad(likfun,p0,dat=ddat,times=t0,N=N)  ## finite differences
## matches!

library(bbmle)
parnames(likfun) <- names(p0)

m1 <- mle2(likfun,start=p0,gr=gradlikfun,data=list(times=t0,N=N,dat=ddat),
           vecpar=TRUE)

plot(t0[-1],ddat)
lines(t0[-1],calc_mean(coef(m1),times=t0,N=N))

pp1 <- profile(m1,which="lambda")

m0 <- mle2(likfun,start=p0,data=list(times=t0,N=N,dat=ddat),
           vecpar=TRUE)

pp0 <- profile(m0,which="lambda")
par(mfrow=c(1,2))
plot(pp1,show.points=TRUE)
plot(pp0,show.points=TRUE)
