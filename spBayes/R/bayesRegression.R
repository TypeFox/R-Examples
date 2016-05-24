##set.seed(1234)

#library(MASS)
##library(mvtnorm)
##library(fields)
##library(magic)
#library(coda)

##takes args like MASS's mvrnorm
rmvn <- function(n, mu=0, V = matrix(1))
{
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  matrix(rnorm(n*p), ncol=p) %*% D + rep(mu,rep(n,p))
}

##final.sampler <- function (dataframe, totalsample, burnin, thinning)
##{
##select <- seq(from=burnin+1, to=totalsample, by=thinning)
##temp <- dataframe[select,]
##temp
##}

##gibbs.sample <- function(file, burnin, thinning)
##{
##dataframe <- read.table(file, header=F)
##select <- seq(from=burnin+1, to = nrow(dataframe), by=thinning)
##output <- dataframe[select,]
##output <- as.matrix(output)
##output
##}

##param.quantiles <- function(gibbsout)
##{
##
##  if (is.data.frame(gibbsout)) {
##    p <- ncol(gibbsout)
##    temp <- matrix(rep(0,times=p*3), ncol=3, byrow=T)
##
##	for(i in 1:p)
##	{
##	temp[i,] <- quantile(gibbsout[,i], c(0.50, 0.025, 0.975))
##      }
##  }
##  else {
##    p <- 1
##    temp <- quantile(gibbsout, c(0.50,0.025,0.975))
##  }
##  
##  temp
##}

bayesNormalReference <- function(sample.mean, sample.var, sample.size, NITER) {
  n = sample.size
  y.bar = sample.mean
  s.sq = sample.var
  sigmasq.post <- 1/rgamma(NITER, (n-1)/2, (n-1)*s.sq/2)
  mu.post <- rep(0,times=NITER)
  for (i in 1:NITER) {
    mu.post[i] = rnorm(1, y.bar, sqrt(sigmasq.post[i]/n))
  }
  out = cbind(mu.post, sigmasq.post)
  out = as.data.frame(out)
  names(out) = c("mu.posterior","sigmasq.posterior")
  return(out)
}


##bayesLMRef <- function(lm.obj, NITER) {
##  lm.obj.data = model.frame(lm.obj)
##  lm.obj.cov = vcov(lm.obj)
##  ##D.inv = diag(diag(solve(lm.objcov)))
##  ##unit.corr = sqrt(D.inv)%*%unit.cov%*%sqrt(D.inv)
##  coeff.mean = coefficients(lm.obj)
##
##  ##unit.coeff.post = rmvn(1000, mu=unit.coeff.mean, Sigma=unit.cov)
##  ##unit.coeff.post.qnt = param.quantiles(unit.coeff.post)
##
##  sigmasq = (summary(lm.obj)$sigma)^2
##
##  N = nrow(lm.obj.data)
##  p = ncol(lm.obj.cov)
##  ## Simulate from marginal posterior [sigma^2 | y]
##  sigmasq.post <- 1/rgamma(NITER, (N-p)/2, (N-p)*sigmasq/2)
##  
##  ## Simulate from conditional posterior [beta | sigma^2, y]
##  coeff.post <- matrix(0, nrow=NITER,ncol=p)
##  for (i in 1:NITER) {
##    Sigma.coeff = (sigmasq.post[i]/sigmasq)*lm.obj.cov
##    ##    coeff.post[i,] = rmvn(1, coeff.mean, Sigma.coeff)
##    coeff.post[i,] = rmvnorm(1, coeff.mean, Sigma.coeff)
##  }
##
##  posterior.samples = as.data.frame(cbind(coeff.post, sigmasq.post))
##  
##  Parameter.summary = data.frame(param.quantiles(posterior.samples),row.names=c(names(coefficients(lm.obj)),"sigma^2")) 
##  names(Parameter.summary) =  c("50%","2.5%","97.5%")
##  Parameter.summary
##}

## ##same as above, but with some cosmetic changes
## bayesLMRef <- function(lm.obj, n.samples, ...) {

##   ####################################################
##   ##Check for unused args
##   ####################################################
##   formal.args <- names(formals(sys.function(sys.parent())))
##   elip.args <- names(list(...))
##   for(i in elip.args){
##     if(! i %in% formal.args)
##       warning("'",i, "' is not an argument")
##   }

##   if(missing(lm.obj)){stop("error: lm.obj must be specified")}
##   if(class(lm.obj) != "lm"){stop("error: lm.obj must be of class lm")}
##   if(missing(n.samples)){stop("error: n.samples must be specified")} 
  
##   lm.obj.data <- model.frame(lm.obj)
##   lm.obj.cov <- vcov(lm.obj)
##   ##D.inv = diag(diag(solve(lm.objcov)))
##   ##unit.corr = sqrt(D.inv)%*%unit.cov%*%sqrt(D.inv)
##   coeff.mean <- coefficients(lm.obj)

##   ##unit.coeff.post = rmvn(1000, mu=unit.coeff.mean, Sigma=unit.cov)
##   ##unit.coeff.post.qnt = param.quantiles(unit.coeff.post)

##   sigmasq <- (summary(lm.obj)$sigma)^2

##   N <- nrow(lm.obj.data)
##   p <- ncol(lm.obj.cov)
##   ## Simulate from marginal posterior [sigma^2 | y]
##   sigmasq.post <- 1/rgamma(n.samples, (N-p)/2, (N-p)*sigmasq/2)
  
##   ## Simulate from conditional posterior [beta | sigma^2, y]
##   coeff.post <- matrix(0, nrow=n.samples,ncol=p)
##   for (i in 1:n.samples) {
##     Sigma.coeff <- (sigmasq.post[i]/sigmasq)*lm.obj.cov
##     coeff.post[i,] <- rmvn(1, coeff.mean, Sigma.coeff)
##   }

##   Parameter.post <- cbind(coeff.post, sigmasq.post)
##   colnames(Parameter.post) <- c(names(coefficients(lm.obj)),"sigma.sq")
##   ##mcmc(Parameter.post)
##   ##Parameter.summary <- data.frame(param.quantiles(Parameter.post),row.names=c(names(coefficients(lm.obj)),"sigma.sq")) 
##   ##names(Parameter.summary) <- c("50%","2.5%","97.5%")
##   ##Parameter.summary


##   ##
##   ##return
##   ##
##   out <- list()
##   out$p.samples <- mcmc(Parameter.post)
##   out$X <- as.matrix(model.matrix(lm.obj))
##   out$Y <- as.matrix(model.extract(model.frame(lm.obj), "response"))
##   out$n.samples <- n.samples
##   class(out) <- "bayesLMRef"
##   out
## }
 

normalIGammaSampler <- function(n.samples, mu, V, a, b){

  sigmasq <-  1/rgamma(n.samples, a, b)

  p <- length(mu)
  beta <- matrix(0, nrow=n.samples,ncol=p)

  for (i in 1:n.samples) {
##    beta[i,] <- rmvn(1, mu, sigmasq[i]*V)
##    beta[i,] <- rmvnorm(1, mu, sigmasq[i]*V)
    beta[i,] <- rmvn(1, mu, sigmasq[i]*V)
  }

  Output <- as.data.frame(cbind(beta, sigmasq))
  
  Output
}

bayesLMConjugate <- function(formula, data = parent.frame(), n.samples, beta.prior.mean, beta.prior.precision,
                               prior.shape, prior.rate, ...) {

  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }

  if(missing(n.samples)){stop("error: n.samples must be specified")} 
  if(missing(beta.prior.mean)){stop("error: beta.prior.mean must be specified")}
  if(missing(beta.prior.precision)){stop("error: beta.prior.precision must be specified")}
  if(missing(prior.shape)){stop("error: prior.shape must be specified")}
  if(missing(prior.rate)){stop("error: prior.rate must be specified")}

  
  if(missing(formula)){stop("error: formula must be specified")} 
  if(class(formula) == "formula"){
    
    holder <- parseFormula(formula, data)
    Y <- holder[[1]]
    X <- as.matrix(holder[[2]])
    x.names <- holder[[3]]

  }else{
    stop("error: formula is misspecified")
  }

  p <- ncol(X)
  n <- nrow(X)


  ##check for dim problems
  if(length(beta.prior.mean) != p)
    stop(paste("error: beta.prior.mean must be of length ",p, sep=""))

  beta.prior.precision <- as.matrix(beta.prior.precision)
  if(nrow(beta.prior.precision) != ncol(beta.prior.precision))
    stop("error: beta.prior.precision must be square")

  if(nrow(beta.prior.precision) != p)
    stop(paste("error: beta.prior.precision must be of dimension ",p, sep=""))    
  
    
  V.beta.inv <- beta.prior.precision
  mu <- beta.prior.mean
  
  tXX <- crossprod(X)
  
  posterior.precision <- V.beta.inv + tXX
  posterior.variance <- chol2inv(chol(posterior.precision))
  posterior.mean <- posterior.variance%*%(V.beta.inv%*%mu + t(X)%*%Y)

  posterior.shape <- prior.shape + n/2
  posterior.rate <- prior.rate + 0.5*(t(mu)%*%V.beta.inv%*%mu + sum(Y*Y) - t(posterior.mean)%*%posterior.precision%*%posterior.mean)

  posterior.samples <- as.matrix(normalIGammaSampler(n.samples, posterior.mean, posterior.variance, posterior.shape, posterior.rate))
  colnames(posterior.samples) <- c(x.names,"sigma.sq")
  ##mcmc(posterior.samples)
  
  
  #Parameter.summary = data.frame(param.quantiles(posterior.samples),row.names=c(names(coefficients(lm.obj)),"sigma.sq")) 
  #names(Parameter.summary) =  c("50%","2.5%","97.5%")
  #Parameter.summary

  ##
  ##return
  ##
  out <- list()
  out$p.samples <- mcmc(posterior.samples)
  out$X <- as.matrix(X)
  out$Y <- as.matrix(Y)
  out$n.samples <- n.samples
  class(out) <- "bayesLMConjugate"
  out
}

##bayes.geostat.exp.unmarg <- function (lm.obj, n.samples, beta.prior.mean, beta.prior.precision, coords, phi, alpha,
##                                      nugget.prior.shape, nugget.prior.rate) {
##  
##  lm.obj.data = model.frame(lm.obj)
##  
##  n <- nrow(lm.obj.data)
##  
##  Y <- lm.obj.data[,1]
##  X <- lm.obj.data[,2:ncol(lm.obj.data)]
##  X <- as.matrix(cbind(rep(1, times=nrow(lm.obj.data)), X))
##
##  p <- ncol(X)
##
##  tildeX <- cbind(X, diag(n)) 
##  ttildeXtildeX <- crossprod(tildeX)
##  
##  D <- as.matrix(dist(coords))
##  R.exp <- exp(-phi*D)
##  spatial.prior.precision <- chol2inv(chol(R.exp))
##  
##  V.theta.inv <- (1/alpha)* adiag(beta.prior.precision, spatial.prior.precision)
##  mu <- c(beta.prior.mean,rep(0,times=n))
##  
##  posterior.precision <- V.theta.inv + ttildeXtildeX
##  posterior.variance <- chol2inv(chol(posterior.precision))
##  posterior.mean <- posterior.variance%*%(V.theta.inv%*%mu + t(tildeX)%*%Y)
##
##  nugget.posterior.shape <- prior.shape + n/2
##  nugget.posterior.rate <- prior.rate + 0.5*(t(mu)%*%V.theta.inv%*%mu + sum(Y*Y) - t(posterior.mean)%*%posterior.precision%*%posterior.mean)
##
##  posterior.samples <- normalIGammaSampler(n.samples, posterior.mean, posterior.variance, nugget.posterior.shape, nugget.posterior.rate)
##
##  beta.posterior.samples <- posterior.samples[,1:p]
##  spatial.effects.posterior.samples <- posterior.samples[,(p+1):(p+n)]
##  nugget.posterior.samples <- posterior.samples[,(n+p+1)]
##  spatial.variance.posterior.samples <- alpha*nugget.posterior.samples
## 
##  list("beta.p.samples"=mcmc(beta.posterior.samples),
##       "spatial.effects.p.samples"=mcmc(spatial.effects.posterior.samples),
##       "nugget.p.samples"=mcmc(nugget.posterior.samples),
##       "spatial.variance.p.samples"=mcmc(spatial.variance.posterior.samples))
##
##}

## NOTE: The following function is a *faster* marginalized version of the above function.
## WARNING: The alpha here is tau^2/sigma^2 (noise-to-signal) and is the reciprocal of what
## was used in the above function.

bayesGeostatExact <- function (formula, data = parent.frame(), n.samples, beta.prior.mean, beta.prior.precision,
                                 coords, cov.model="exponential", phi, nu, alpha, sigma.sq.prior.shape, sigma.sq.prior.rate,
                                 sp.effects=TRUE, verbose=TRUE, ...) {

  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }

  if(missing(n.samples)){stop("error: n.samples must be specified")} 
  if(missing(beta.prior.mean)){stop("error: beta.prior.mean must be specified")}
  if(missing(beta.prior.precision)){stop("error: beta.prior.precision must be specified")}
  if(missing(coords)){stop("error: coords must be specified")}
  if(missing(phi)){stop("error: phi must be specified")}
  if(missing(alpha)){stop("error: alpha must be specified")}
  if(missing(sigma.sq.prior.shape)){stop("error: sigma.sq.prior.shape must be specified")}
  if(missing(sigma.sq.prior.rate)){stop("error: sigma.sq.prior.rate must be specified")}

  if(!cov.model%in%c("gaussian","exponential","matern","spherical"))
    {stop("error: specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")}

  if(cov.model == "matern")
    if(missing(nu)){stop("error: nu must be specified")}
  
  if(missing(formula)){stop("error: formula must be specified")} 
  if(class(formula) == "formula"){
    
    holder <- parseFormula(formula, data)
    Y <- holder[[1]]
    X <- as.matrix(holder[[2]])
    x.names <- holder[[3]]

  }else{
    stop("error: formula is misspecified")
  }

  p <- ncol(X)
  n <- nrow(X)

  ##check for dim problems
  if(length(beta.prior.mean) != p)
    stop(paste("error: beta.prior.mean must be of length ",p, sep=""))

  beta.prior.precision <- as.matrix(beta.prior.precision)
  if(nrow(beta.prior.precision) != ncol(beta.prior.precision))
    stop("error: beta.prior.precision must be square")

  if(nrow(beta.prior.precision) != p)
    stop(paste("error: beta.prior.precision must be of dimension ",p, sep=""))    

  ##
  ##show summary of stuff
  ##
  if(verbose){
    cat("-------------------------------------------------\n")
    cat("\tGeneral model description\n")
    cat("-------------------------------------------------\n")
    cat(paste("Model fit with ",n," observations.\n", sep=""))
    cat(paste("Number of covariates ",p," (including intercept if specified).\n", sep=""))
    cat(paste("Using the ",cov.model," spatial correlation model.\n\n", sep=""))
    cat("-------------------------------------------------\n")
    cat("\t\tSampling\n")
    cat("-------------------------------------------------\n")
  }


  
  D <- as.matrix(dist(coords))

  if(cov.model == "exponential"){
    R <- exp(-phi*D)
  }else if(cov.model == "matern"){
    R <- (D*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=D*phi, nu=nu)
    diag(R) <- 1
  }else if(cov.model == "gaussian"){
    R <- exp(-1*((phi*D)^2))
  }else if(cov.model == "spherical"){
    R <- D
    R[TRUE] <- 1
    R[D > 0 & D < 1/phi] <- 1-1.5*phi*D[D > 0 & D <= 1/phi]+0.5*((phi*D[D > 0 & D <= 1/phi])^3)
    R[D >= 1/phi] <- 0   
  }else{
    stop("error: specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")
  }

  V.y <- R + alpha*diag(nrow(R))   
  V.y.sqrt <- chol(V.y)  
  V.y.inv <- chol2inv(V.y.sqrt)
  ttildeXtildeX <- t(X)%*%V.y.inv%*%X
  mu <- beta.prior.mean
  
  posterior.precision <- beta.prior.precision + ttildeXtildeX
  posterior.variance <- chol2inv(chol(posterior.precision))
  posterior.mean <- posterior.variance%*%(beta.prior.precision%*%mu + t(X)%*%V.y.inv%*%Y)

  sigma.sq.posterior.shape <- sigma.sq.prior.shape + n/2
  sigma.sq.posterior.rate <- sigma.sq.prior.rate + 0.5*(t(mu)%*%beta.prior.precision%*%mu + t(Y)%*%V.y.inv%*%Y - t(posterior.mean)%*%posterior.precision%*%posterior.mean)

  posterior.samples <- as.matrix(normalIGammaSampler(n.samples, posterior.mean, posterior.variance, sigma.sq.posterior.shape, sigma.sq.posterior.rate))

  beta.posterior.samples <- posterior.samples[,1:p]
  sigma.sq.posterior.samples <- posterior.samples[,(p+1)]
  tau.sq.posterior.samples <- alpha*sigma.sq.posterior.samples  

  posterior.samples <- cbind(beta.posterior.samples, sigma.sq.posterior.samples, tau.sq.posterior.samples)
  colnames(posterior.samples) <- c(x.names,"sigma.sq", "tau.sq")

  cat(paste("Sampled: ",n.samples," of ",n.samples,", ",100,"%\n", sep=""))
  
  if(sp.effects){
    cat("-------------------------------------------------\n")
    cat("\tRecovering spatial effects\n")
    cat("-------------------------------------------------\n")

    w <- matrix(0, n, n.samples)

    R.inv <- chol2inv(chol(R))
    
    V.sp <- chol2inv(chol(R.inv + (1/alpha)*diag(nrow(R))))
    resid.posterior <- matrix(rep(Y, times=n.samples), nrow=length(Y), ncol=n.samples) -  X%*%t(beta.posterior.samples) 
    sp.posterior.mean <- (1/alpha)*t(V.sp%*%resid.posterior)
    
    V.sp.root <- t(chol(V.sp)) ## chol returns "upper-triangular"; so t(); 
    
    
    for (s in 1:n.samples) {
      
      ## Using rmvnorm is slow - it calculates the matrix square-root each time
      ## sp.effects[s,] <- rmvnorm(1, sp..posterior.mean[s,], sigma.sq.posterior.samples[s]*V.sp)
      
      ## Instead use the pre-computed V.sp.root
      z <- rnorm(nrow(V.sp), 0, 1)
      w[,s] <- sp.posterior.mean[s,] + sqrt(sigma.sq.posterior.samples[s])*V.sp.root%*%z
      
    }
    
##    w <- matrix(0, n, n.samples)
##    
##    R.eigen <- eigen(R)
##    R.vals.inv <- 1/R.eigen$values
##
##    R.vecs <- R.eigen$vectors
##
##    sigma.sq <- posterior.samples[,"sigma.sq"]
##    tau.sq <- posterior.samples[,"tau.sq"]
##    beta <- as.matrix(posterior.samples[,1:p])
##
##    R.vects.t <- t(R.vecs)
##    
##    for(s in 1:n.samples){
##      
##      S.w <- R.vecs%*%diag(1/(1/sigma.sq[s]*R.vals.inv+1/tau.sq[s]))%*%R.vects.t
##      
##      S.mu <- S.w%*%(Y - X%*%as.matrix(beta[s,]))/tau.sq[s]
##      
##      S.w.sq <- R.vecs%*%diag(sqrt(1/(1/sigma.sq[s]*R.vals.inv+1/tau.sq[s])))
##
##      w[,s] <- S.w.sq%*%as.matrix(rnorm(n))+S.mu
##    }
    
  } ##end sp.effects


  ##
  ##return
  ##
  out <- list()

  out$X <- X
  out$n <- n
  out$p <- p
  out$Y <- Y
  out$n.samples <- n.samples
  out$beta.prior.mean <- beta.prior.mean
  out$beta.prior.precision <- beta.prior.precision
  out$coords <- coords
  out$cov.model <- cov.model
  out$phi <- phi
  out$alpha <- alpha
  out$sigma.sq.prior.shape <- sigma.sq.prior.shape
  out$sigma.sq.prior.rate <- sigma.sq.prior.rate
  out$alpha <- alpha
  out$verbose <-verbose
 
  if(cov.model == "matern")
    out$nu <- nu
  
  out$p.samples <- mcmc(posterior.samples)
  
  if(sp.effects)
    out$sp.effects <- w

  class(out) <- "bayesGeostatExact"
  
  out
}

