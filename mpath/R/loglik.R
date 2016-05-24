### log-likelihood function, except for gaussian family
logLik.glmreg <- function(object, newx, y, weights, na.action=na.pass, ...){
if(missing(newx) || missing(y)) return(object$twologlik/2)
if(!is.null(object$terms)){
 mf <- model.frame(delete.response(object$terms), newx, na.action = na.action, xlev = object$xlevels)
 newx <- model.matrix(delete.response(object$terms), mf, contrasts = object$contrasts)
 newx <- newx[,-1] ### remove the intercept
}
family <- object$family
nlambda <- object$nlambda
beta <- object$beta
b0 <- object$b0
famtype <- switch(family,
      "gaussian"=1,
      "binomial"=2,
      "poisson"=3,
      "negbin"=4)
 nm <- dim(newx)
 nobs <- nm[1]
 nvars <- nm[2]
 if(is.matrix(beta)){
 if(dim(beta)[1] != nvars)
 stop("newx has a wrong dimension, possibly the 1st column for the intecept\n")
}
if(missing(weights)) weights <- rep(1, nobs)
w <- weights
 res <- .Fortran("pred",
 n=as.integer(nobs),
 m=as.integer(nvars),
 nlambda = as.integer(nlambda),
 x=as.double(newx),
 b=as.double(beta),
 a0=as.double(b0),
 family=as.integer(famtype),
 eta=as.double(matrix(0, nobs, nlambda)),
 mu=as.double(matrix(0, nobs, nlambda)),
 package="mpath")
 mu <- matrix(res$mu, ncol=nlambda) 
 th <- object$theta
 if(nlambda > 1 && length(object$theta) ==1)
 th <- rep(th, nlambda)
# res <- llfun(y, mu, w, th=th, family, nlambda)
 res <- rep(NA, nlambda)
 for(k in 1:nlambda) 
 res[k] <- .Fortran("loglikFor",
  n=as.integer(length(y)), 
  y=as.double(y), 
  mu=as.double(mu[,k]), 
  theta=as.double(th[k]), 
  w=as.double(w), 
  family=as.integer(famtype),
  ll=as.double(0),
  package="mpath")$ll
 return(res)
}

llfun <- function(y, mu, w, th, family, nlambda){
 ll <- rep(NA, nlambda)
switch(family,
  "gaussian"={
    for(k in 1:nlambda)
    ll[k] <- -sum(w*(y-mu[,k])^2)
},
  "negbin"={
    for(k in 1:nlambda)
    ll[k] <- sum(w*(lgamma(th[k] + y) - lgamma(th[k]) - lgamma(y + 1) + th[k] * log(th[k]) +
           y * log(mu[,k] + (y == 0)) - (th[k] + y) * log(th[k] + mu[,k])))
},
  "binomial"={
    for(k in 1:nlambda){
     tmp <- rep(0, length(y))
    for(i in 1:length(y))
     if(is.na(mu[i,k])) cat("mu[i,k]=", mu[i,k], "\n")
     if(mu[i,k]>0 && mu[i,k] <1)
      tmp[i] <- w[i]*(y[i]*log(mu[i,k]) + (1-y[i])*log(1-mu[i,k]))
      ll[k] <- sum(tmp)
}
},
  "poisson"={
    for(k in 1:nlambda)
    ll[k] <- sum(w*(-mu[,k] + y*log(mu[,k])-lgamma(y+1)))
}
)
return(ll)
}

AIC.glmreg <- function(object, ..., k)
object$aic
BIC.glmreg <- function(object, ...)
object$bic
AIC.zipath <- function(object, ..., k)
object$aic
BIC.zipath <- function(object, ...)
object$bic

logLik.zipath <- function(object, newdata, y, weights, na.action=na.pass, 
                   link = c("logit", "probit", "cloglog", "cauchit", "log"),
...){
if(missing(newdata) || missing(y)) return(object$loglik)
linkstr <- match.arg(link)
  linkobj <- make.link(linkstr)
  linkinv <- linkobj$linkinv
Y <- y
## set up likelihood
  ziPoisson <- function(parms) {
  kx <- NCOL(X)
  kz <- NCOL(Z)
    Y0 <- Y <= 0
    Y1 <- Y > 0
    offsetx <- offsetz <- 0
    ## count mean
    mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
    ## binary mean
    phi <- as.vector(linkinv(Z %*% parms[(kx+1):(kx+kz)] + offsetz))

    ## log-likelihood for y = 0 and y >= 1
    loglik0 <- log( phi + exp( log(1-phi) - mu ) ) ## -mu = dpois(0, lambda = mu, log = TRUE)
    loglik1 <- log(1-phi) + dpois(Y, lambda = mu, log = TRUE)

    ## collect and return
    loglik <- sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] * loglik1[Y1])
    loglik
  }

  ziNegBin <- function(parms) {
  kx <- NCOL(X)
  kz <- NCOL(Z)
    Y0 <- Y <= 0
    Y1 <- Y > 0
    offsetx <- offsetz <- 0
    ## count mean
    mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
    ## binary mean
    phi <- as.vector(linkinv(Z %*% parms[(kx+1):(kx+kz)] + offsetz))
    ## negbin size
    theta <- exp(parms[(kx+kz)+1])

    ## log-likelihood for y = 0 and y >= 1
    loglik0 <- log( phi + exp( log(1-phi) + suppressWarnings(dnbinom(0, size = theta, mu = mu, log = TRUE)) ) )
    loglik1 <- log(1-phi) + suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))

    ## collect and return
    loglik <- sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] * loglik1[Y1])
    loglik
  }
  ziGeom <- function(parms) 
  ziNegBin(c(parms, 0))

family <- object$family
nlambda <- object$nlambda
  loglikfun <- switch(family,
                      "poisson" = ziPoisson,
                      "geometric" = ziGeom,
                      "negbin" = ziNegBin)

ll <- rep(NA, nlambda)
 mf <- model.frame(delete.response(object$terms$full), newdata, na.action = na.action, xlev = object$levels)
    X <- model.matrix(delete.response(object$terms$count), mf, contrasts = object$contrasts$count)
    Z <- model.matrix(delete.response(object$terms$zero),  mf, contrasts = object$contrasts$zero)
if(missing(weights)) weights <- rep(1, length(y))
w <- weights
if(family == "negbin")
parms <- rbind(object$coefficients$count, object$coefficients$zero, log(object$theta))
else parms <- rbind(object$coefficients$count, object$coefficients$zero)
 ll <- rep(NA, nlambda)
for(k in 1:nlambda)
ll[k] <- loglikfun(parms[,k])
ll
}
