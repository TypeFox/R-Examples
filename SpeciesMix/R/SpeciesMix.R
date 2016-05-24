# This is package SpeciesMix 

"additive.logistic" <-
function (x,inv=FALSE) 
{
  if(inv){
    x <- log(x/x[length(x)])
    return(x)
  }

  x.t <- exp(x)
  x.t <- x.t/(1+sum(x.t))
  x.t[length(x.t)+1] <- 1-sum(x.t)
  return(x.t)
}


"apply.glm" <-
function (i,form,datsp,tau,n) 
{
  dat.tau <- rep(tau[,i],each=n)
  x <- model.matrix(as.formula(form),data=datsp)
  y <- datsp$obs
  ##dat.tau <- get("dat.tau")
  ##datsp <- get("datsp")
  ##f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau[,i],family="binomial",x=T,y=T)
  
  f.mix <- glm.fit(x,y,dat.tau,family=binomial())
  return(list(coef=f.mix$coef))
 
 
  ##f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau,family="binomial",x=T,y=T)
  ##list(residuals=f.mix$residuals,fitted=f.mix$fitted,linear.predictors=f.mix$linear.predictors,coef=f.mix$coef)

}


"apply.glm.gaussian" <-
function (i,form,datsp,tau,n) 
{
  dat.tau <- rep(tau[,i],each=n)
  x <- model.matrix(as.formula(form),data=datsp)
  y <- datsp$obs
  environment(form) <- environment()
  ##dat.tau <- get("dat.tau")
  ##datsp <- get("datsp")
  ##f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau[,i],family="binomial",x=T,y=T)

 
    ##f.mix <- glm.nb(as.formula(form),data=datsp,weights=dat.tau)
  f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau,family="gaussian")
  ##f.mix <- glm.nbinom(as.formula(form),data=datsp,weights=dat.tau)
  sp.int <- rep(f.mix$coef[1],dim(tau)[1])
    ##return(list(coef=f.mix$coef[-1],theta=f.mix$theta))
  return(list(coef=f.mix$coef[-1],theta=sqrt(f.mix$deviance/f.mix$df.residual),sp.intercept=sp.int))
  ##return(list(coef=f.mix$coef,theta=f.mix$theta))
 
  ##f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau,family="binomial",x=T,y=T)
  ##list(residuals=f.mix$residuals,fitted=f.mix$fitted,linear.predictors=f.mix$linear.predictors,coef=f.mix$coef)

}


"apply.glm.nbinom" <-
function (i,form,datsp,tau,n) 
{
  dat.tau <- rep(tau[,i],each=n)
  x <- model.matrix(as.formula(form),data=datsp)
  y <- datsp$obs
  environment(form) <- environment()
  ##dat.tau <- get("dat.tau")
  ##datsp <- get("datsp")
  ##f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau[,i],family="binomial",x=T,y=T)

 
    ##f.mix <- glm.nb(as.formula(form),data=datsp,weights=dat.tau)
  #f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau,family="poisson")
  f.mix <- glm.nbinom(as.formula(form),data=datsp,weights=dat.tau)
  sp.int <- rep(f.mix$coef[1],dim(tau)[1])
    ##return(list(coef=f.mix$coef[-1],theta=f.mix$theta))
  return(list(coef=f.mix$coef[-1],theta=f.mix$theta,sp.intercept=sp.int))
  ##return(list(coef=f.mix$coef,theta=f.mix$theta))
 
  ##f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau,family="binomial",x=T,y=T)
  ##list(residuals=f.mix$residuals,fitted=f.mix$fitted,linear.predictors=f.mix$linear.predictors,coef=f.mix$coef)

}


"apply.glm.tweedie" <-
function (i,form,datsp,tau,n) 
{
  dat.tau <- rep(tau[,i],each=n)
  x <- model.matrix(as.formula(form),data=datsp)
  y <- datsp$obs
  environment(form) <- environment()
   ##f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau[,i],family="binomial",x=T,y=T)

 
  f.mix <- tglm(as.formula(form),data=datsp,wts=dat.tau,vcov=FALSE,residuals=FALSE,trace=0)
  sp.int <- rep(f.mix$coef[1],dim(tau)[1])
 
  return(list(coef=f.mix$coef[c(-1,-(length(f.mix$coef)-1),-(length(f.mix$coef)))],phi=f.mix$coef["phi"],p=f.mix$coef["p"],sp.intercept=sp.int))
 
}


"artificial.data" <-
function (formula, data, theta, S, dist = "bernoulli") 
{
    X <- model.matrix(formula, data)
    out <- matrix(0, dim(X)[1], S)
    k <- dim(theta)[1]
    sp.int <- rep(0, S)
    group <- rep(0, S)
    for (s in 1:S) {
        g <- ceiling(runif(1) * k)
        if (dist == "bernoulli") {
            lgtp <- X %*% theta[g, ]
            p <- exp(lgtp)/(1 + exp(lgtp))
            out[, s] <- rbinom(dim(X)[1], 1, p)
        }
        if (dist == "negbin") {
            tmp <- rep(1e+05, dim(X)[1])
            while (max(tmp, na.rm = T) > 5000 | sum(tmp) < 100) {
                theta[g, 1] <- runif(1, -15, 5)
                sp.int[s] <- theta[g, 1]
                lgtp <- X %*% theta[g, ]
                p <- exp(lgtp)
                tmp <- rnbinom(dim(X)[1], mu = p, size = 1)
            }
            out[, s] <- tmp
        }
        if (dist == "tweedie") {
            tmp <- rep(6e+05, dim(X)[1])
            while (max(tmp, na.rm = T) > 5e+05 | sum(tmp) < 100) {
                theta[g, 1] <- runif(1, -15, 5)
                theta[g, 1] <- runif(1, 1, 5)
                sp.int[s] <- theta[g, 1]
                lgtp <- X %*% theta[g, ]
                p <- exp(lgtp)
                tmp <- rTweedie(dim(X)[1], mu = p, phi = 2, p = 1.6)
            }
            out[, s] <- tmp
        }
        if (dist == "gaussian") {
            tmp <- rep(1e+05, dim(X)[1])
            while (max(tmp, na.rm = T) > 50000 | sum(tmp) < 100) {
                theta[g, 1] <- runif(1, 100, 500)
                lgtp <- X %*% theta[g, ]
                p <- (lgtp)
                tmp <- rnorm(dim(X)[1], mean = p, sd = 1)
            }
            out[, s] <- tmp
        }
        group[s] <- g
    }
    pi <- tapply(group, group, length)/S
    if (dist == "negbin") 
        return(list(pa = out, group = group, pi = pi, sp.int = sp.int))
    if (dist == "tweedie") 
        return(list(pa = out, group = group, pi = pi, sp.int = sp.int))
    list(pa = out, group = group, pi = pi)
}


"clusterSelect" <-
function (sp.form,sp.data,covar.data,G=1:10,em.prefit=TRUE, em.steps=4 ,em.refit=3,est.var=FALSE,trace=TRUE) 
{
  my.fun <- function(g,form,sp.data,covar.data){
    cat("Fitting group",g,"\n")
    try(SpeciesMix(form,sp.data,covar.data,g,em.prefit=em.prefit,em.steps=em.steps,em.refit=em.refit,est.var=est.var,trace=trace))
  }
  

#  if(mc){ out <- mclapply(G,my.fun,form,dat,mc.preschedule = FALSE, mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = set.cores)} else
  { out <- lapply(G,my.fun,sp.form,sp.data,covar.data)}

  aic <- rep(0,length(G))
  bic <- rep(0,length(G))
  fm <- list()
  for(i in 1:length(G))
    if(!is.atomic(out[[i]])){
      aic[i] <- out[[i]]$aic
      bic[i] <- out[[i]]$bic
      fm[[i]] <- list(logl=out[[i]]$logl,coef=out[[i]]$coef,tau=out[[i]]$tau,pi=out[[i]]$pi,covar=out[[i]]$covar)
    }
  return(list(aic=aic,bic=bic,fm=fm))
 
}


"create.starting.values" <-
function (S,G,n,form,datsp) 
{
  ##apply.glm <- function(i,form,datsp,tau,n,dat.tau){
  ##  dat.tau <- rep(tau[,i],each=n)
  ##  dat.tau <- get("dat.tau")
  ##  datsp <- get("datsp")
    ##f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau[,i],family="binomial",x=T,y=T)
  ##  f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau,family="binomial",x=T,y=T)
    ##list(residuals=f.mix$residuals,fitted=f.mix$fitted,linear.predictors=f.mix$linear.predictors,coef=f.mix$coef)
   ## list(coef=f.mix$coef)
  ##}
  
  environment(form) <- environment()
  tau <- matrix(runif(S*G),S,G)
  tau <- (tau/rowSums(tau))
  fmM <- list()
  for(i in 1:G){
      ##dat.tau[,i] <- rep(tau[,i],each=n)
      pi[i] <- sum(tau[,i])/S
    }
  ##dat.tau <- rep(0,dim(datsp)[1])

  fmM <- lapply(1:G,apply.glm,form,datsp,tau,n)
  first.fit <- list(x=model.matrix(as.formula(form),data=datsp),y=datsp$obs,formula=form)
 
  ##return(list(pi=pi,fmM=fmM,tau=tau,dat.tau=dat.tau,first.fit=first.fit))
  return(list(pi=pi,fmM=fmM,tau=tau,first.fit=first.fit))
}


"create.starting.values.gaussian" <-
function (S,G,n,form,datsp,mc=FALSE,set.cores=2) 
{
  ##apply.glm <- function(i,form,datsp,tau,n,dat.tau){
  ##  dat.tau <- rep(tau[,i],each=n)
  ##  dat.tau <- get("dat.tau")
  ##  datsp <- get("datsp")
    ##f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau[,i],family="binomial",x=T,y=T)
  ##  f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau,family="binomial",x=T,y=T)
    ##list(residuals=f.mix$residuals,fitted=f.mix$fitted,linear.predictors=f.mix$linear.predictors,coef=f.mix$coef)
   ## list(coef=f.mix$coef)
  ##}
  
  environment(form) <- environment()
  tau <- matrix(runif(S*G),S,G)
  tau <- (tau/rowSums(tau))
  fmM <- list()
  for(i in 1:G){
      ##dat.tau[,i] <- rep(tau[,i],each=n)
      pi[i] <- sum(tau[,i])/S
    }
  ##dat.tau <- rep(0,dim(datsp)[1])
 
  fmM <- lapply(1:G,apply.glm.gaussian,form,datsp,tau,n)
  first.fit <- list(x=model.matrix(as.formula(form),data=datsp)[,-1],y=datsp$obs,formula=form)
 
  ##return(list(pi=pi,fmM=fmM,tau=tau,dat.tau=dat.tau,first.fit=first.fit))
  return(list(pi=pi,fmM=fmM,tau=tau,first.fit=first.fit))
}


"create.starting.values.nbinom" <-
function (S,G,n,form,datsp,mc=FALSE,set.cores=2) 
{
  ##apply.glm <- function(i,form,datsp,tau,n,dat.tau){
  ##  dat.tau <- rep(tau[,i],each=n)
  ##  dat.tau <- get("dat.tau")
  ##  datsp <- get("datsp")
    ##f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau[,i],family="binomial",x=T,y=T)
  ##  f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau,family="binomial",x=T,y=T)
    ##list(residuals=f.mix$residuals,fitted=f.mix$fitted,linear.predictors=f.mix$linear.predictors,coef=f.mix$coef)
   ## list(coef=f.mix$coef)
  ##}
  
  environment(form) <- environment()
  tau <- matrix(runif(S*G),S,G)
  tau <- (tau/rowSums(tau))
  fmM <- list()
  for(i in 1:G){
      ##dat.tau[,i] <- rep(tau[,i],each=n)
      pi[i] <- sum(tau[,i])/S
    }
  ##dat.tau <- rep(0,dim(datsp)[1])
 
  fmM <- lapply(1:G,apply.glm.nbinom,form,datsp,tau,n)
  offset <- model.frame(as.formula(form),data=datsp)
  offset <- model.offset(offset)
  if(is.null(offset)) offset <- rep(0,length(datsp$obs))
  first.fit <- list(x=model.matrix(as.formula(form),data=datsp)[,-1],y=datsp$obs,offset=offset,formula=form)
 
  ##return(list(pi=pi,fmM=fmM,tau=tau,dat.tau=dat.tau,first.fit=first.fit))
  return(list(pi=pi,fmM=fmM,tau=tau,first.fit=first.fit))
}


"create.starting.values.nbinom.kmeans" <-
function (S, G, n, form, datsp, tol = 0.1) 
{
    MM <- model.matrix(form, datsp)
    offset <- model.frame(form,datsp)
    offset <- model.offset(offset)
    if(is.null(offset)) offset <- rep(0,nrow(MM))
    
    first.fit <- list(x = model.matrix(as.formula(form), data = datsp)[, 
        -1], y = datsp$obs, offset=offset, formula = form)
    if (tol < 0 || tol >= 1) 
        stop("Minimum Prevalence % must be between 0 and 1")
    sp.name <- 1:S
    sp <- rep(sp.name, each = n)
    starting.fitem <- list(intercepts = rep(0, S), alpha = rep(0, 
        S))
    all.betas <- matrix(0, nrow = S, ncol = ncol(MM))
    colnames(all.betas) <- colnames(MM)
    for (j in 1:S) {
        fit <- glm.nbinom(form, datsp[((j - 1) * n):(j * n - 
            1), ], est.var = FALSE)
        starting.fitem$sp.intercepts[j] <- fit$coef[1]
        all.betas[j, ] <- fit$coef
        starting.fitem$theta[j] = fit$theta
    }
    all.betas[, 1] <- 0
    cat("Clustering...\n")
    fmmvnorm <- kmeans(x = all.betas, centers = G, iter.max = 100, 
        nstart = 50)
    starting.fitem$coef <- fmmvnorm$centers
    fmM <- list()
    for (i in 1:G) {
        B <- matrix(rep(fmmvnorm$centers[i, ], nrow(datsp)), 
            nrow(datsp), ncol(fmmvnorm$centers), byrow = T)
        B[, 1] <- rep(starting.fitem$sp.intercepts, each = n)
        fitted <- exp(rowSums(MM * B)+offset)
        fmM[[i]] <- list(coef = fmmvnorm$centers[i, 2:ncol(fmmvnorm$centers)], 
            theta = mean(starting.fitem$theta[fmmvnorm$cluster == 
                i]), sp.intercept = starting.fitem$sp.intercepts, 
            fitted = fitted)
    }
    tau <- matrix(0, S, G)
    pi <- rep(1/G, G)
    pi <- runif(G, 0.2, 0.8)
    pi <- pi/sum(pi)
    est.tau <- lapply(1:S, estimate.pi.nbinom, sp, sp.name, datsp, 
        fmM, pi, G, first.fit)
    max.newTau <- 0.8
    alpha <- (1 - max.newTau * G)/(max.newTau * (2 - G) - 1)
    for (j in 1:S) {
        newTau <- (2 * alpha * est.tau[[j]]$tau - alpha + 1)/(2 * 
            alpha - alpha * G + G)
        tau[j, ] <- newTau
    }
    for (i in 1:G) {
        pi[i] <- sum(tau[, i])/S
    }
    return(list(pi = pi, fmM = fmM, tau = tau, first.fit = first.fit))
}


"create.starting.values.tweedie" <-
function (S,G,n,form,datsp) 
{
  ##apply.glm <- function(i,form,datsp,tau,n,dat.tau){
  ##  dat.tau <- rep(tau[,i],each=n)
  ##  dat.tau <- get("dat.tau")
  ##  datsp <- get("datsp")
    ##f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau[,i],family="binomial",x=T,y=T)
  ##  f.mix <- glm(as.formula(form),data=datsp,weights=dat.tau,family="binomial",x=T,y=T)
    ##list(residuals=f.mix$residuals,fitted=f.mix$fitted,linear.predictors=f.mix$linear.predictors,coef=f.mix$coef)
   ## list(coef=f.mix$coef)
  ##}
  
  environment(form) <- environment()
  tau <- matrix(runif(S*G),S,G)
  tau <- (tau/rowSums(tau))
  fmM <- list()
  for(i in 1:G){
      ##dat.tau[,i] <- rep(tau[,i],each=n)
      pi[i] <- sum(tau[,i])/S
    }
  ##dat.tau <- rep(0,dim(datsp)[1])
  offset <- model.frame(form,datsp)
  offset <- model.offset(offset)
  if(is.null(offset)) offset <- rep(0,length(datsp$obs))
  fmM <- lapply(1:G,apply.glm.tweedie,form,datsp,tau,n)
  first.fit <- list(x=model.matrix(as.formula(form),data=datsp)[,-1],y=datsp$obs,formula=form)
 
  ##return(list(pi=pi,fmM=fmM,tau=tau,dat.tau=dat.tau,first.fit=first.fit))
  return(list(pi=pi,fmM=fmM,tau=tau,first.fit=first.fit))
}


"create.starting.values.tweedie.kmeans" <-
function (S, G, n, form, datsp, tol = 0.1) 
{
    MM <- model.matrix(form, datsp)
     offset <- model.frame(form,datsp)
    offset <- model.offset(offset)
    if(is.null(offset)) offset <- rep(0,nrow(MM))
    first.fit <- list(x = model.matrix(as.formula(form), data = datsp)[, 
        -1], y = datsp$obs, offset=offset,formula = form)
    if (tol < 0 || tol >= 1) 
        stop("Minimum Prevalence % must be between 0 and 1")
    sp.name <- 1:S
    sp <- rep(sp.name, each = n)
    starting.fitem <- list(intercepts = rep(0, S), alpha = rep(0, 
        S))
    all.betas <- matrix(0, nrow = S, ncol = ncol(MM))
    colnames(all.betas) <- colnames(MM)
    for (j in 1:S) {
        fit <- tglm(form, datsp[((j - 1) * n):(j * n - 1), ], 
            vcov = FALSE, residuals = FALSE, trace = 0, p = 1.6)
        starting.fitem$sp.intercepts[j] <- fit$coef[1]
        all.betas[j, ] <- fit$coef[-length(fit$coef)]
        starting.fitem$phi[j] = fit$coef["phi"]
    }
    all.betas[, 1] <- 0
    cat("Clustering...\n")
    fmmvnorm <- kmeans(x = all.betas, centers = G, iter.max = 100, 
        nstart = 50)
    starting.fitem$coef <- fmmvnorm$centers
    fmM <- list()
    for (i in 1:G) {
        B <- matrix(rep(fmmvnorm$centers[i, ], nrow(datsp)), 
            nrow(datsp), ncol(fmmvnorm$centers), byrow = T)
        B[, 1] <- rep(starting.fitem$sp.intercepts, each = n)
        fitted <- exp(rowSums(MM * B)+offset)
        fmM[[i]] <- list(coef = fmmvnorm$centers[i, 2:ncol(fmmvnorm$centers)], 
            phi = mean(starting.fitem$phi[fmmvnorm$cluster == 
                i]), p = 1.6, sp.intercept = starting.fitem$sp.intercepts, 
            fitted = fitted)
    }
    tau <- matrix(0, S, G)
    pi <- rep(1/G, G)
    pi <- runif(G, 0.2, 0.8)
    pi <- pi/sum(pi)
    est.tau <- lapply(1:S, estimate.pi.tweedie, sp, sp.name, 
        datsp, fmM, pi, G, first.fit)
    max.newTau <- 0.8
    alpha <- (1 - max.newTau * G)/(max.newTau * (2 - G) - 1)
    for (j in 1:S) {
        newTau <- (2 * alpha * est.tau[[j]]$tau - alpha + 1)/(2 * 
            alpha - alpha * G + G)
        tau[j, ] <- newTau
    }
    for (i in 1:G) {
        pi[i] <- sum(tau[, i])/S
    }
    return(list(pi = pi, fmM = fmM, tau = tau, first.fit = first.fit))
}


"distr.binom" <-
function( p){
 nobs <- length( p)
  new.dist <- old.dist <- rep( 0, nobs+1)
  old.dist[1] <- 1-p[1]
  old.dist[2] <- p[1]
  for( ii in 2:nobs){
    new.dist[1] <- old.dist[1]*(1-p[ii])
    for( jj in 2:ii)
      new.dist[jj] <- old.dist[jj-1]*p[ii] + old.dist[jj]*(1-p[ii])
    new.dist[ii+1] <- old.dist[ii]*p[ii]
    old.dist <- new.dist
  }
  return( new.dist)
}


"dPoisGam" <-
function ( y, lambda, mu.Z, alpha, LOG=TRUE) 
{
#function to calculate Random sum (Tweedie) densities.
#y is the value of the r.v.  Can be a vector
#mu.N is the mean of the Poisson summing r.v. Can be a vector of length(y)
#mu.Z is the mean of the Gamma rv Can be a vector of length(y)
#alpha is the `other' parameter of the gamma distribution s.t. var = ( mu.Z^2)/alpha Can be a vector of length(y)
#If mu.N, mu.Z or alpha are scalare but y isn't then they will be used for all y. If lengths mis-match then error
#LOG=TRUE gives the density on the log scale
#do.checks=TRUE checks the input vectors for compatability and gives errors / changes them as appropriate.
#do.checks=FALSE doesn't check and relies on the user to have things right. If not right then catastrophic failure may occur.

#  if( any( is.null( c( y, mu.N, mu.Z, alpha)))){
#    print( "Error: null input values -- please check.  Null values are:")
#    tmp <- double( is.null( c( y, mu.N, mu.Z, alpha)))
#    names( tmp) <- c( "y", "mu.N","mu.Z","alpha")
#    print( tmp)
#    print( "Exitting")
#    return()
#  }
  mu.N <- lambda
  if( !all( is.element( c( length( mu.N), length( mu.Z), length( alpha)), c( length( y), 1)))){
    print( "Error: length of parameter vectors does not match length of random variable vector")
    print( "Exitting")
    return()
  }

  if( length( mu.N) != length( y))
    mu.N <- rep( mu.N, length( y))
  if( length( mu.Z) != length( y))
    mu.Z <- rep( mu.Z, length( y))
  if( length( alpha) != length( y))
    alpha <- rep( alpha, length( y))

  res <- .Call( "dTweedie", as.numeric( y), as.numeric( mu.N), as.numeric( mu.Z), as.numeric( alpha), as.integer( LOG),PACKAGE="SpeciesMix")

  return( res)

}


"dPoisGamDerivs" <-
function ( y=NULL, lambda=NULL, mu.Z=NULL, alpha=NULL, do.checks=TRUE) 
{
#function to calculate Random sum (Tweedie) densities.
#y is the value of the r.v.  Can be a vector
#mu.N is the mean of the Poisson summing r.v. Can be a vector of length(y)
#mu.Z is the mean of the Gamma rv Can be a vector of length(y)
#alpha is the `other' parameter of the gamma distribution s.t. var = ( mu.Z^2)/alpha Can be a vector of length(y)
#If mu.N, mu.Z or alpha are scalare but y isn't then they will be used for all y. If lengths mis-match then error
#LOG=TRUE gives the density on the log scale
#do.checks=TRUE checks the input vectors for compatability and gives errors / changes them as appropriate.
#do.checks=FALSE doesn't check and relies on the user to have things right. If not right then catastrophic failure may occur.

  mu.N <- lambda
  if( do.checks){
    if( any( is.null( c( y, mu.N, mu.Z, alpha)))){
      print( "Error: null input values -- please check.  Null values are:")
      tmp <- double( is.null( c( y, mu.N, mu.Z, alpha)))
      names( tmp) <- c( "y", "mu.N","mu.Z","alpha")
      print( tmp)
      print( "Exitting")
      return()
    }

    if( !all( is.element( c( length( mu.N), length( mu.Z), length( alpha)), c( length( y), 1)))){
      print( "Error: length of parameter vectors does not match length of random variable vector")
      print( "Exitting")
    }

    if( length( mu.N) != length( y))
      mu.N <- rep( mu.N, length( y))
    if( length( mu.Z) != length( y))
      mu.Z <- rep( mu.Z, length( y))
    if( length( alpha) != length( y))
      alpha <- rep( alpha, length( y))
  }

  res <- .Call( "dTweedieDeriv", as.numeric( y), as.numeric( mu.N), as.numeric( mu.Z), as.numeric( alpha),PACKAGE="SpeciesMix")
  colnames( res) <- c("lambda","mu.Z","alpha")
  return( res)

}


"dTweedie" <-
function ( y, mu, phi, p, LOG=TRUE) 
{
  lambda <- ( mu^( 2-p)) / ( phi*(2-p))
  alpha <- ( 2-p) / ( p-1)
  tau <- phi*(p-1)*mu^(p-1)
  mu.Z <- alpha * tau

  dens <- dPoisGam( y, lambda, mu.Z, alpha, LOG)
  
  return( dens)

}


"estimate.pi" <-
function (j,sp,spname,datsp,fmM,pi,G,first.fit) 
{
   
  tmp.like <- rep(0,G)
  tau <- rep(0,G)
  sel.sp <- which(sp==spname[j])
  
  link.fun <- make.link("logit")
  for(i in 1:G) {
    if(length(fmM[[i]]$coef)==1){lpre <- link.fun$linkinv(first.fit$x[sel.sp,]*fmM[[i]]$coef)
                               }else{    lpre <- link.fun$linkinv(first.fit$x[sel.sp,]%*%fmM[[i]]$coef)}
    ##lpre <- link.fun$linkinv(first.fit$x[sel.sp,]%*%fmM[[i]]$coef)
    ## tmp.like[i] <- sum(dbinom(datsp$obs[sel.sp],1,fmM[[i]]$fitted[sel.sp],log=T))
    obs <- datsp$obs[sel.sp]
    lpre[obs==0]<- 1- lpre[obs==0]
    tmp.like[i] <- sum(log(lpre))
  }
  eps <- max(tmp.like)
  sum.like <- (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
  for(i in 1:G) {
    tau[i] <- exp((log(pi[i]) + tmp.like[i]) - sum.like )
  }

  return(list(tau=tau,sum.like=sum.like))
}


"estimate.pi.gaussian" <-
function (j,sp,spname,datsp,fmM,pi,G,first.fit) 
{
   
  tmp.like <- rep(0,G)
  tau <- rep(0,G)
  sel.sp <- which(sp==spname[j])
    
 
  for(i in 1:G) {
    ##lpre <- dnorm(first.fit$y[sel.sp],mean=fmM[[i]]$fitted[sel.sp],sd=fmM[[i]]$theta,log=TRUE)
    lpre <- dnorm(first.fit$y[sel.sp],mean=fmM[[i]]$fitted[sel.sp],sd=1,log=TRUE)
    ##lpre <- dnbinom(first.fit$y[sel.sp],mu=exp(cbind(1,first.fit$x[sel.sp,])%*%c(fmM[[i]]$sp.intercept[j],fmM[[i]]$coef)),size=fmM[[i]]$theta,log=TRUE)
     ##lpre <- dpois(first.fit$y[sel.sp],exp(cbind(1,first.fit$x[sel.sp,])%*%c(fmM[[i]]$sp.intercept[j],fmM[[i]]$coef)),log=TRUE)
    #lpre <- dpois(first.fit$y[sel.sp],exp(first.fit$x[sel.sp,]%*%fmM[[i]]$coef),log=TRUE)
    
    tmp.like[i] <- sum(lpre)
  }
  
    
  eps <- max(tmp.like)
  sum.like <- (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
  for(i in 1:G) {
    tau[i] <- exp((log(pi[i]) + tmp.like[i]) - sum.like )
  }
  return(list(tau=tau,sum.like=sum.like))
}


"estimate.pi.nbinom" <-
function (j,sp,spname,datsp,fmM,pi,G,first.fit) 
{
   
  tmp.like <- rep(0,G)
  tau <- rep(0,G)
  sel.sp <- which(sp==spname[j])
    
 
  for(i in 1:G) {
    lpre <- dnbinom(first.fit$y[sel.sp],mu=fmM[[i]]$fitted[sel.sp],size=fmM[[i]]$theta,log=TRUE)
    ##lpre <- dnbinom(first.fit$y[sel.sp],mu=exp(cbind(1,first.fit$x[sel.sp,])%*%c(fmM[[i]]$sp.intercept[j],fmM[[i]]$coef)),size=fmM[[i]]$theta,log=TRUE)
     ##lpre <- dpois(first.fit$y[sel.sp],exp(cbind(1,first.fit$x[sel.sp,])%*%c(fmM[[i]]$sp.intercept[j],fmM[[i]]$coef)),log=TRUE)
    #lpre <- dpois(first.fit$y[sel.sp],exp(first.fit$x[sel.sp,]%*%fmM[[i]]$coef),log=TRUE)
    
    tmp.like[i] <- sum(lpre)
  }
  
    
  eps <- max(tmp.like)
  sum.like <- (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
  for(i in 1:G) {
    tau[i] <- exp((log(pi[i]) + tmp.like[i]) - sum.like )
  }
  return(list(tau=tau,sum.like=sum.like))
}


"estimate.pi.tweedie" <-
function (j, sp, spname, datsp, fmM, pi, G, first.fit) 
{
    tmp.like <- rep(0, G)
    tau <- rep(0, G)
    sel.sp <- which(sp == spname[j])
    for (i in 1:G) {
        lpre <- dTweedie(first.fit$y[sel.sp], mu = fmM[[i]]$fitted[sel.sp], 
            phi = fmM[[i]]$phi, p = fmM[[i]]$p, LOG = TRUE)
        tmp.like[i] <- sum(lpre)
    }
    eps <- max(tmp.like)
    sum.like <- (log(sum(pi * exp((tmp.like) - (eps)))) + (eps))
    for (i in 1:G) {
        tau[i] <- exp((log(pi[i]) + tmp.like[i]) - sum.like)
    }
    return(list(tau = tau, sum.like = sum.like))
}


"fitMix" <-
function (form,datsp,sp,G=2,ite.max=500,trace=TRUE,full.model=FALSE,r1=FALSE) 
{
## dat2 has colums obs,sp
  ##
  temp.warn <- getOption( "warn")
  options( warn=-1)

   sp.name <- unique(sp)
  S <- length(unique(sp))
  n <- length(which(sp==sp.name[1]))

  cat("Fitting Group",G,"\n")
  if(trace) cat("Iteration | LogL \n")
  
  ##dat.tau <- data.frame(matrix(0,dim(datsp)[1],G))
  dat.tau <- 0
  ##dat <- data.frame(datsp,dat.tau)
  pi <- rep(0,G)
  ite <- 1
  logL <- -99999999
  old.logL <- -88888888

  ## set up initial GLM
  t1 <- create.starting.values(S,G,n,form,datsp)
  pi <- t1$pi;fmM <- t1$fmM;tau <- t1$tau;dat.tau <- t1$dat.tau;first.fit <- t1$first.fit
 
 
  while(abs(logL-old.logL) > 0.0001 & ite<=ite.max){
    old.logL <- logL
    
    for(i in 1:G){
      ##dat.tau[,i] <- rep(tau[,i],each=n)
      pi[i] <- sum(tau[,i])/S
    }
    
    if(any(pi==0)) { ## occasionally with complicated models the random starting values result in a pi[i]==0; so restart with new random starts
      cat("pi has gone to zero - restarting fitting \n")
      t1 <- create.starting.values(S,G,n,form,datsp)
      pi <- t1$pi;fmM <- t1$fmM;tau <- t1$tau;dat.tau <- t1$dat.tau;first.fit <- t1$first.fit
      ite <- 1
    }

    fmM <- lapply(1:G,weighted.glm,first.fit,tau,n,fmM,sp) 


    logL <- 0 
    tmp.like <- matrix(0,S,G)

    est.tau <- lapply(1:S,estimate.pi,sp,sp.name,datsp,fmM,pi,G,first.fit)
           
    for(j in 1:S){
      if(is.atomic(est.tau[[j]])){ print (est.tau[[j]])} else
      {
        tau[j,] <- est.tau[[j]]$tau
        logL <- logL+est.tau[[j]]$sum.like
      }
    }

    if(trace) cat(ite,"     | ",logL,"\n")
     ite <- ite+1
  }
  fm.out <- data.frame(matrix(0,G,length(fmM[[1]]$coef)))
  names(fm.out) <- names(fmM[[1]]$coef)
  tau <- data.frame(tau)
  names(tau) <- paste("grp.",1:G,sep="")
  EN <- -sum(unlist(tau)*log(unlist(tau)))
  d <- length(unlist(fm.out)) + length(tau)-1
  for(i in 1:G) {
    fm.out[i,] <- fmM[[i]]$coef
    ##dat.tau[,i] <- rep(tau[,i],each=n)
  }

  names(pi) <- paste("G",1:G,sep=".")
  t.pi <- additive.logistic(pi,TRUE)
  parms <- c(t.pi[1:(G-1)],unlist(fm.out))
  logL.full <- logL
  logL <- logLmix(parms,first.fit,G,S,sp,sp.name)

  options(warn=temp.warn)
  if(full.model)  return(list(logl=logL,aic = -2*logL + 2*d,tau=round(tau,4),pi=pi,bic=-2*logL + log(S)*d,ICL= -2*logL + log(S)*d +2*EN,coef=fm.out,fmM=fmM,model.tau=dat.tau,covar=0,aic.full= -2*logL.full + 2*d,bic.full= -2*logL.full + log(S)*d))
  
  return(list(logl=logL,aic = -2*logL + 2*d,tau=round(tau,4),pi=pi,bic=-2*logL + log(S)*d,ICL= -2*logL + log(S)*d +2*EN,coef=fm.out,covar=0,aic.full= -2*logL.full + 2*d,bic.full= -2*logL.full + log(S)*d))
  
}


"fitmix.cpp" <-
function (form,datsp,sp,G=2,pars=NA,trace=TRUE,calc.hes=FALSE,r1) 
{
  if(!is.numeric(sp)){
    sp <- as.integer(factor(sp))
  }
  sp.name <- unique(sp)
  S <- length(unique(sp))
  n <- length(which(sp==sp.name[1]))

  X <- model.matrix(form, data = datsp[sp==sp.name[1],])
  ##X <- model.frame(form, data = datsp)
 # X <- X[sp==sp.name[1],]
  
  #y <- model.response(form, data = datsp)
  y <- datsp$obs
  if(is.na(pars[1])) {
    ##pars <- rep(0.01,G-1+(ncol(X)*G))
    pars <- runif(G-1+(ncol(X)*G),-2,2)
    pars[1:(G-1)] <- runif(G-1)
  }
  offset<-rep(0,n)
##  hes <- rep(0,length(pars)^2)
  gradient <- rep(0,length(pars))
  tau <- matrix(0,S,G) ##must leave this in as defines S & G
  model.type<- as.integer(1)
  loglike <- try(.Call("SpeciesMix",pars,y,X,sp,tau,gradient,offset,model.type,PACKAGE="SpeciesMix"))
  
  calc.deriv <- function(p){
    gradient <- rep(0,length(pars))
    ll <- .Call("Calculate_Gradient",p,y,X,sp,tau,gradient,offset,as.integer(1),PACKAGE="SpeciesMix")
    return(gradient)
  }
  hes <- 0
  if(calc.hes){
    hes <- nd2(pars,calc.deriv)
    dim(hes) <- rep(length(pars),2)
    dim(hes) <- rep(length(pars),2)
    rownames(hes) <- colnames(hes) <- c(paste("G.",1:(G-1),sep=""),paste("G",1:G,rep(colnames(X),each=G),sep="."))
  }
  if(!is.numeric(loglike)) loglike <- 0
  pi <- pars[1:(G-1)]
  coef <- pars[ (G):length(pars)]

  r.logl <- logLmix(pars,list(y=y,x=model.matrix(form, data = datsp)),G,S,sp,sp.name,out.tau=TRUE)
  pi <- additive.logistic(pi)
  names(pi) <- paste("G.",1:G,sep="")
  coef <- matrix(coef,G,ncol(X))
  rownames(coef) <- paste("G.",1:G,sep="")
  colnames(coef) <- colnames(X)

  AIC <- 2*loglike + 2*length(pars)
  BIC <- 2*loglike + log(S)*length(pars)
  ##list(logl=loglike,r.logl=r.logl$logl,pi=pi,coef=coef,tau=round(exp(r.logl$tau),4),aic=AIC,bic=BIC,hessian=hes,gradient=gradient)
  list(logl=loglike,pi=pi,coef=coef,tau=round(exp(r.logl$tau),4),aic=AIC,bic=BIC,hessian=hes,gradient=gradient)
  
}


"fitmix.gaussian" <-
function (form,datsp,sp,G=2,ite.max=500,trace=TRUE,full.model=FALSE) 
{
## fitting abundance data with a negative binomial distribution
  ## dat2 has colums obs,sp
  ##
  temp.warn <- getOption( "warn")
  options( warn=-1)

  sp.name <- unique(sp)
  S <- length(unique(sp))
  n <- length(which(sp==sp.name[1]))

  cat("Fitting Group",G,"\n")
  if(trace) cat("Iteration | LogL \n")
  
  ##dat.tau <- data.frame(matrix(0,dim(datsp)[1],G))
  dat.tau <- 0
  ##dat <- data.frame(datsp,dat.tau)
  pi <- rep(0,G)
  ite <- 1
  logL <- -99999999
  old.logL <- -88888888

  ## set up initial GLM
  t1 <- create.starting.values.gaussian(S,G,n,form,datsp)
  pi <- t1$pi;fmM <- t1$fmM;tau <- t1$tau;dat.tau <- t1$dat.tau;first.fit <- t1$first.fit
 
 
  while(abs(logL-old.logL) > 0.0001 & ite<=ite.max){
    old.logL <- logL
    
    for(i in 1:G){
      ##dat.tau[,i] <- rep(tau[,i],each=n)
      pi[i] <- sum(tau[,i])/S
    }
    
    if(any(pi==0)) { ## occasionally with complicated models the random starting values result in a pi[i]==0; so restart with new random starts
      cat("pi has gone to zero - restarting fitting \n")
      t1 <- create.starting.values.gaussian(S,G,n,form,datsp)
      pi <- t1$pi;fmM <- t1$fmM;tau <- t1$tau;dat.tau <- t1$dat.tau;first.fit <- t1$first.fit
      ite <- 1
    }

     fmM <- lapply(1:G,weighted.glm.gaussian,first.fit,tau,n,fmM,sp) 


    logL <- 0 
    tmp.like <- matrix(0,S,G)

    est.tau <- lapply(1:S,estimate.pi.gaussian,sp,sp.name,datsp,fmM,pi,G,first.fit)
           
    for(j in 1:S){
      if(is.atomic(est.tau[[j]])){ print (est.tau[[j]])} else
      {
        tau[j,] <- est.tau[[j]]$tau
        logL <- logL+est.tau[[j]]$sum.like
      }
    }

    if(trace) cat(ite,"     | ",logL,"\n")
     ite <- ite+1
  }
  fm.out <- data.frame(matrix(0,G,length(fmM[[1]]$coef)))
  int.out <- data.frame(matrix(0,S,G))
  names(int.out) <- paste("grp.",1:G,sep="")
  names(fm.out) <- names(fmM[[1]]$coef)
  tau <- data.frame(tau)
  names(tau) <- paste("grp.",1:G,sep="")
  EN <- -sum(unlist(tau)*log(unlist(tau)))
  d <- length(unlist(fm.out)) + length(tau)-1
  fm.theta <- rep(0,G)
  for(i in 1:G) {
    fm.out[i,] <- fmM[[i]]$coef
    int.out[,i] <- fmM[[i]]$sp.intercept
    ##dat.tau[,i] <- rep(tau[,i],each=n)
    fm.theta[i] <- fmM[[i]]$theta
  }

  names(pi) <- paste("G",1:G,sep=".")
  t.pi <- additive.logistic(pi,TRUE)
  parms <- c(t.pi[1:(G-1)],fm.theta,unlist(fm.out),unlist(int.out))
  logL.full <- logL
  logL <- logLmix.gaussian(parms,first.fit,G,S,sp,sp.name)

  options(warn=temp.warn)
  if(full.model)  return(list(logl=logL,aic = -2*logL + 2*d,tau=round(tau,4),pi=pi,bic=-2*logL + log(S)*d,ICL= -2*logL + log(S)*d +2*EN,coef=fm.out,sp.intercept=int.out,theta=fm.theta,fmM=fmM,model.tau=dat.tau,covar=0,aic.full= -2*logL.full + 2*d,bic.full= -2*logL.full + log(S)*d))
  
  return(list(logl=logL,aic = -2*logL + 2*d,tau=round(tau,4),pi=pi,bic=-2*logL + log(S)*d,ICL= -2*logL + log(S)*d +2*EN,coef=fm.out,sp.intercept=int.out,theta=fm.theta,covar=0,aic.full= -2*logL.full + 2*d,bic.full= -2*logL.full + log(S)*d))
  
}


"fitmix.nbinom" <-
function (form, datsp, sp, G = 2, ite.max = 500, trace = TRUE, 
    full.model = FALSE) 
{
    temp.warn <- getOption("warn")
    options(warn = -1)
    sp.name <- unique(sp)
    S <- length(unique(sp))
    n <- length(which(sp == sp.name[1]))
    cat("Fitting Group", G, "\n")
    if (trace) 
        cat("Iteration | LogL \n")
    dat.tau <- 0
    pi <- rep(0, G)
    ite <- 1
    logL <- -99999999
    old.logL <- -88888888
    if (ite != 1) 
        t1 <- create.starting.values.nbinom(S, G, n, form, datsp)
    t1 <- create.starting.values.nbinom.kmeans(S, G, n, form, 
        datsp)
    pi <- t1$pi
    fmM <- t1$fmM
    tau <- t1$tau
    dat.tau <- t1$dat.tau
    first.fit <- t1$first.fit
    while (abs(logL - old.logL) > 1e-04 & ite <= ite.max) {
        old.logL <- logL
        for (i in 1:G) {
            pi[i] <- sum(tau[, i])/S
        }
        if (any(pi == 0)) {
            cat("pi has gone to zero - restarting fitting \n")
            t1 <- create.starting.values.nbinom(S, G, n, form, 
                datsp)
            pi <- t1$pi
            fmM <- t1$fmM
            tau <- t1$tau
            dat.tau <- t1$dat.tau
            first.fit <- t1$first.fit
            ite <- 1
        }
        fmM <- lapply(1:G, weighted.glm.nbinom, first.fit, tau, 
            n, fmM, sp)
        for (j in 1:S) {
            tmp <- rep(0, G)
            for (g in 1:G) tmp[g] <- fmM[[g]]$sp.intercept[j]
            tmp <- sum(tmp * tau[j, ])
            for (g in 1:G) fmM[[g]]$sp.intercept[j] <- tmp
        }
        logL <- 0
        tmp.like <- matrix(0, S, G)
        est.tau <- lapply(1:S, estimate.pi.nbinom, sp, sp.name, 
            datsp, fmM, pi, G, first.fit)
        for (j in 1:S) {
            if (is.atomic(est.tau[[j]])) {
                print(est.tau[[j]])
            }
            else {
                tau[j, ] <- est.tau[[j]]$tau
                logL <- logL + est.tau[[j]]$sum.like
            }
        }
        if (trace) 
            cat(ite, "     | ", logL, "\n")
        ite <- ite + 1
    }
    fm.out <- data.frame(matrix(0, G, length(fmM[[1]]$coef)))
    int.out <- rep(0, S)
    names(fm.out) <- names(fmM[[1]]$coef)
    tau <- data.frame(tau)
    names(tau) <- paste("grp.", 1:G, sep = "")
    EN <- -sum(unlist(tau) * log(unlist(tau)))
    d <- length(unlist(fm.out)) + length(tau) - 1
    fm.theta <- rep(0, G)
    int.out <- fmM[[1]]$sp.intercept
    for (i in 1:G) {
        fm.out[i, ] <- fmM[[i]]$coef
    }
    names(pi) <- paste("G", 1:G, sep = ".")
    t.pi <- additive.logistic(pi, TRUE)
    parms <- c(t.pi[1:(G - 1)], unlist(fm.out), int.out, rep(1, 
        S))
    logL.full <- logL
    logL <- logLmix.nbinom(parms, first.fit, G, S, sp, sp.name)
    options(warn = temp.warn)
    if (full.model) 
        return(list(logl = logL, aic = -2 * logL + 2 * d, tau = round(tau, 
            4), pi = pi, bic = -2 * logL + log(S) * d, ICL = -2 * 
            logL + log(S) * d + 2 * EN, coef = fm.out, sp.intercept = int.out, 
            theta = fm.theta, fmM = fmM, model.tau = dat.tau, 
            covar = 0, aic.full = -2 * logL.full + 2 * d, bic.full = -2 * 
                logL.full + log(S) * d, pars = parms))
    return(list(logl = logL, aic = -2 * logL + 2 * d, tau = round(tau, 
        4), pi = pi, bic = -2 * logL + log(S) * d, ICL = -2 * 
        logL + log(S) * d + 2 * EN, coef = fm.out, sp.intercept = int.out, 
        theta = fm.theta, covar = 0, aic.full = -2 * logL.full + 
            2 * d, bic.full = -2 * logL.full + log(S) * d, pars = parms))
}


"fitmix.nbinom.cpp" <-
function (form,datsp,sp,G=2,pars=NA,trace=TRUE,calc.hes=FALSE) 
{
  if(!is.numeric(sp)){
    sp <- as.integer(factor(sp))
  }
  sp.name <- unique(sp)
  S <- length(unique(sp))
  n <- length(which(sp==sp.name[1]))

  X <- model.matrix(form, data = datsp[sp==sp.name[1],])
##  offset <- model.offset(form, data = datsp[sp==sp.name[1],])
  offset <- model.frame(form, data = datsp[sp==sp.name[1],])
  offset <- model.offset(offset)
  if(is.null(offset)) offset <-  rep(0,n)
  ##offset <-  rep(0,n)
  ##X <- model.frame(form, data = datsp)
 # X <- X[sp==sp.name[1],]
  
  #y <- model.response(form, data = datsp)
  y <- datsp$obs
  if(is.na(pars[1])) {
    ##pars <- rep(0.01,G-1+(ncol(X)*G))
    sp.int <- rep(0.5,S)
    sp.dispersion <- rep(1,S)
    fm <- matrix(runif(ncol(X)*G,-1,1),G,ncol(X))
    pars <- c(runif(G-1),unlist(fm),sp.int,sp.dispersion)
    
  }

##  hes <- rep(0,length(pars)^2)
  gradient <- rep(0,length(pars))
  tau <- matrix(0,S,G) ##must leave this in as defines S & G
  
  loglike <- try(.Call("SpeciesMix",pars,y,X,sp,tau,gradient,offset,as.integer(2),PACKAGE="SpeciesMix"))
  
  calc.deriv <- function(p){
    gradient <- rep(0,length(pars))
    ll <- .Call("Calculate_Gradient",p,y,X,sp,tau,gradient,offset,as.integer(2),PACKAGE="SpeciesMix")
    return(gradient)
  }
  r.deriv <- function(p){ logLmix.nbinom(p,list(y=y,x=model.matrix(form, data = datsp)),G,S,sp,sp.name,out.tau=FALSE)}
  #r.grad <- nd2(pars,r.deriv)
  ##print(r.grad)
  hes <- 0
  covar <- 0
  if(calc.hes){
    hes <- nd2(pars,calc.deriv)
    dim(hes) <- rep(length(pars),2)
    dim(hes) <- rep(length(pars),2)
    covar <- try(solve(hes))
    #rownames(hes) <- colnames(hes) <- c(paste("G.",1:(G-1),sep=""),paste("G",1:G,rep(colnames(X),each=G),sep="."))
  }
  if(!is.numeric(loglike)) loglike <- 0
  pi <- pars[1:(G-1)]
  #coef <- pars[ (G):length(pars)]
  coef <- pars[-1*(1:((G-1)))]  ## remove pi
  sp.int <- coef[(length(coef)-(2*S-1)):(length(coef)-S)]
  sp.dispersion <- coef[(length(coef)-(S-1)):length(coef)]
  fm <- coef[-1*((length(coef)-(2*S-1)):length(coef))]
  offset <- model.frame(form, data = datsp)
  offset <- model.offset(offset)
  if(is.null(offset)) offset <-  rep(0,nrow(datsp))
  
  r.logl <- logLmix.nbinom(pars,list(y=y,x=model.matrix(form, data = datsp),offset=offset),G,S,sp,sp.name,out.tau=TRUE)
  print(r.logl$logl)
  pi <- additive.logistic(pi)
  names(pi) <- paste("G.",1:G,sep="")
  coef <- matrix(fm,G,ncol(X))
  rownames(coef) <- paste("G.",1:G,sep="")
  colnames(coef) <- colnames(X)

  AIC <- 2*loglike + 2*length(pars)
  BIC <- 2*loglike + log(S)*length(pars)
  ##list(logl=loglike,r.logl=r.logl$logl,pi=pi,coef=coef,tau=round(exp(r.logl$tau),4),aic=AIC,bic=BIC,hessian=hes,gradient=gradient)
  list(logl=loglike,pi=pi,coef=coef,sp.intercept=sp.int,sp.dispersion=sp.dispersion,tau=round(exp(r.logl$tau),4),aic=AIC,bic=BIC,hessian=hes,gradient=gradient,covar=covar)#,r.grad=r.grad)
  
}


"fitmix.tweedie" <-
function (form, datsp, sp, G = 2, ite.max = 500, trace = TRUE, 
    full.model = FALSE) 
{
    temp.warn <- getOption("warn")
    options(warn = -1)
    sp.name <- unique(sp)
    S <- length(unique(sp))
    n <- length(which(sp == sp.name[1]))
    cat("Fitting Group", G, "\n")
    if (trace) 
        cat("Iteration | LogL \n")
    dat.tau <- 0
    pi <- rep(0, G)
    ite <- 1
    logL <- -99999999
    old.logL <- -88888888
    if (ite != 1) 
        t1 <- create.starting.values.tweedie(S, G, n, form, datsp)
    t1 <- create.starting.values.tweedie.kmeans(S, G, n, form, 
        datsp)
    pi <- t1$pi
    fmM <- t1$fmM
    tau <- t1$tau
    dat.tau <- t1$dat.tau
    first.fit <- t1$first.fit
    while (abs(logL - old.logL) > 1e-04 & ite <= ite.max) {
        old.logL <- logL
        for (i in 1:G) {
            pi[i] <- sum(tau[, i])/S
        }
        if (any(pi == 0)) {
            cat("pi has gone to zero - restarting fitting \n")
            t1 <- create.starting.values.tweedie(S, G, n, form, 
                datsp)
            pi <- t1$pi
            fmM <- t1$fmM
            tau <- t1$tau
            dat.tau <- t1$dat.tau
            first.fit <- t1$first.fit
            ite <- 1
        }
        fmM <- lapply(1:G, weighted.glm.tweedie, first.fit, tau, 
            n, fmM, sp)
        for (j in 1:S) {
            tmp <- rep(0, G)
            for (g in 1:G) tmp[g] <- fmM[[g]]$sp.intercept[j]
            tmp <- sum(tmp * tau[j, ])
            for (g in 1:G) fmM[[g]]$sp.intercept[j] <- tmp
        }
        logL <- 0
        tmp.like <- matrix(0, S, G)
        est.tau <- lapply(1:S, estimate.pi.tweedie, sp, sp.name, 
            datsp, fmM, pi, G, first.fit)
        for (j in 1:S) {
            if (is.atomic(est.tau[[j]])) {
                print(est.tau[[j]])
            }
            else {
                tau[j, ] <- est.tau[[j]]$tau
                logL <- logL + est.tau[[j]]$sum.like
            }
        }
        if (trace) 
            cat(ite, "     | ", logL, "\n")
        ite <- ite + 1
    }
    fm.out <- data.frame(matrix(0, G, length(fmM[[1]]$coef)))
    int.out <- fmM[[1]]$sp.intercept
    names(fm.out) <- names(fmM[[1]]$coef)
    tau <- data.frame(tau)
    names(tau) <- paste("grp.", 1:G, sep = "")
    EN <- -sum(unlist(tau) * log(unlist(tau)))
    d <- length(unlist(fm.out)) + length(tau) - 1
    fm.phi <- fm.p <- rep(0, G)
    for (i in 1:G) {
        fm.out[i, ] <- fmM[[i]]$coef
        fm.phi[i] <- fmM[[i]]$phi
        fm.p[i] <- fmM[[i]]$p
    }
    names(pi) <- paste("G", 1:G, sep = ".")
    t.pi <- additive.logistic(pi, TRUE)
    parms <- c(t.pi[1:(G - 1)], unlist(fm.out), int.out, rep(2, 
        S))
    names(parms) <- c(rep("pi", G - 1), rep("coef", length(unlist(fm.out))), 
        rep("int", length(int.out)), rep("phi", S))
    logL.full <- logL
    logL <- logLmix.tweedie(parms, first.fit, G, S, sp, sp.name)
    print(logL)
    options(warn = temp.warn)
    if (full.model) 
        return(list(logl = logL, aic = -2 * logL + 2 * d, tau = round(tau, 
            4), pi = pi, bic = -2 * logL + log(S) * d, ICL = -2 * 
            logL + log(S) * d + 2 * EN, coef = fm.out, sp.intercept = int.out, 
            phi = fm.phi, p = fm.p, fmM = fmM, model.tau = dat.tau, 
            covar = 0, aic.full = -2 * logL.full + 2 * d, bic.full = -2 * 
                logL.full + log(S) * d, pars = parms))
    return(list(logl = logL, aic = -2 * logL + 2 * d, tau = round(tau, 
        4), pi = pi, bic = -2 * logL + log(S) * d, ICL = -2 * 
        logL + log(S) * d + 2 * EN, coef = fm.out, sp.intercept = int.out, 
        phi = fm.phi, p = fm.p, covar = 0, aic.full = -2 * logL.full + 
            2 * d, bic.full = -2 * logL.full + log(S) * d, pars = parms))
}


"fitmix.tweedie.cpp" <-
function (form, datsp, sp, G = 2, pars = NA, trace = TRUE, calc.hes = FALSE) 
{
    if (!is.numeric(sp)) {
        sp <- as.integer(factor(sp))
    }
    sp.name <- unique(sp)
    S <- length(unique(sp))
    n <- length(which(sp == sp.name[1]))
    X <- model.matrix(form, data = datsp[sp == sp.name[1], ])
    offset <- model.frame(form, data = datsp[sp==sp.name[1],])
    offset <- model.offset(offset)
    if(is.null(offset)) offset <-  rep(0,n)

    y <- datsp$obs
    if (is.na(pars[1])) {
        sp.int <- rep(0.5, S)
        sp.phi <- rep(1, S)
        sp.p <- rep(1.6, S)
        fm <- matrix(rep(0.5, ncol(X) * G), G, ncol(X))
        pars <- c(rep(0.5, G - 1), unlist(fm), sp.int, sp.phi)
    }
    pars.og <- pars
    pars.og[1] <- pars.og[1] * 0.99
    gradient <- rep(0, length(pars))
    tau <- matrix(0, S, G)
    loglike <- try(.Call("SpeciesMix", pars, y, X, sp, tau, gradient, 
        offset, as.integer(3),PACKAGE="SpeciesMix"))
    calc.deriv <- function(p) {
        gradient <- rep(0, length(pars))
        ll <- .Call("Calculate_Gradient", p, y, X, sp, tau, gradient, 
            offset, as.integer(3),PACKAGE="SpeciesMix")
        return(gradient)
    }
    hes <- 0
    if (calc.hes) {
        hes <- jacobian(calc.deriv, pars, method = "simple")
        dim(hes) <- rep(length(pars), 2)
        dim(hes) <- rep(length(pars), 2)
    }
    if (!is.numeric(loglike)) 
        loglike <- 0
    pi <- pars[1:(G - 1)]
    coef <- pars[-1 * (1:((G - 1)))]
    sp.int <- coef[(length(coef) - (2 * S - 1)):(length(coef) - 
        S)]
    sp.dispersion <- coef[(length(coef) - (S - 1)):length(coef)]
    sp.phi <- sp.dispersion[1:S]
    sp.p <- rep(1.6, S)
    fm <- coef[-1 * ((length(coef) - (2 * S - 1)):length(coef))]
    offset <- model.frame(form, data = datsp)
    offset <- model.offset(offset)
    if(is.null(offset)) offset <- rep(0,nrow(datsp))
    
    names(gradient) <- names(pars) <- c(rep("pi", G - 1), rep("coef", 
        length(unlist(fm))), rep("int", length(sp.int)), rep("phi", 
        S))
    r.deriv <- function(p) {
        logLmix.tweedie(p, list(y = y, x = model.matrix(form, 
            data = datsp)), G, S, sp, sp.name, out.tau = FALSE)
    }
    r.logl <- logLmix.tweedie(pars, list(y = y, x = model.matrix(form, 
        data = datsp),offset=offset), G, S, sp, sp.name, out.tau = TRUE)
    print(r.logl$logl)
    pi <- additive.logistic(pi)
    names(pi) <- paste("G.", 1:G, sep = "")
    coef <- matrix(fm, G, ncol(X))
    rownames(coef) <- paste("G.", 1:G, sep = "")
    colnames(coef) <- colnames(X)
    AIC <- 2 * loglike + 2 * length(pars)
    BIC <- 2 * loglike + log(S) * length(pars)
    list(logl = loglike, pi = pi, coef = coef, sp.intercept = sp.int, 
        phi = sp.phi, p = sp.p, tau = round(exp(r.logl$tau), 
            4), aic = AIC, bic = BIC, hessian = hes, gradient = gradient)
}


"glm.fit.nbinom" <-
function (x,y,offset=NULL,weights=NULL,mustart=NULL,est.var=FALSE) 
{
  X <- x
  if(is.null(offset)) offset <- 0
  if(is.null(weights)) weights <- rep(1,length(y))
  gradient <- rep(0,ncol(X)+1)
  if(is.null(mustart)){ pars <- gradient+1} else{pars <- mustart}
 fitted.values <- rep(0,length(y))
  logl <- .Call("Neg_Bin",pars,X,y,weights,offset,gradient,fitted.values,PACKAGE="SpeciesMix") 
  vcov <- 0
  se <- rep(0,length(pars))
  if(est.var) {
    calc.deriv <- function(p){
      gradient <- rep(0,length(pars))
      ll <- .Call("Neg_Bin_Gradient",p,X,y,weights,offset,gradient,PACKAGE="SpeciesMix")
      return(gradient)
    }
    hes <- nd2(pars,calc.deriv)
    dim(hes) <- rep(length(pars),2)
    vcov <- try(solve(hes))
    se <- try(sqrt(diag(vcov)))
    colnames(vcov) <- rownames(vcov) <- c("theta",colnames(X))
  }
  names(pars) <- names(se) <- names(gradient) <- c("theta",colnames(X))

  return(list(logl=logl,coef=pars[-1],theta=pars[1],se=se[-1],se.theta=se[1],fitted=fitted.values,gradient=gradient,vcov=vcov))
}


"glm.nbinom" <-
function (form,data,weights=NULL,mustart=NULL,est.var=FALSE) 
{
  X <- model.matrix(form,data)
  t1 <- model.frame(form,data)
  y <- model.response(t1)
  offset <- model.offset(t1)
  if(is.null(offset)) offset <- 0
  if(is.null(weights)) weights <- rep(1,length(y))
  gradient <- rep(0,ncol(X)+1)
  if(is.null(mustart)){ pars <- gradient+1}else{pars <- mustart}
  fitted.values <- rep(0,length(y))

  logl <- .Call("Neg_Bin",pars,X,y,weights,offset,gradient,fitted.values,PACKAGE="SpeciesMix") 
  vcov <- 0
  se <- rep(0,length(pars))
  if(est.var) {
    calc.deriv <- function(p){
      gradient <- rep(0,length(pars))
      ll <- .Call("Neg_Bin_Gradient",p,X,y,weights,offset,gradient,PACKAGE="SpeciesMix")
      return(gradient)
    }
    hes <- nd2(pars,calc.deriv)
    dim(hes) <- rep(length(pars),2)
    vcov <- try(solve(hes))
    se <- try(sqrt(diag(vcov)))
    colnames(vcov) <- rownames(vcov) <- c("theta",colnames(X))
  }
  names(pars) <- names(se) <- names(gradient) <- c("theta",colnames(X))

  return(list(logl=logl,coef=pars[-1],theta=pars[1],se=se[-1],se.theta=se[1],fitted=fitted.values,gradient=gradient,vcov=vcov))
}


"ldTweedie.lp" <-
function ( parms, y, X.p, offsetty, phi, p, wts=rep( 1, length( y))) 
{
  mu <- exp( X.p %*% parms[1:ncol( X.p)] + offsetty)

  if( is.null( phi) & is.null( p)){
    phi <- parms[ncol( X.p) + 1]
    p <- parms[ncol( X.p) + 2]
  }
  if( is.null( phi) & !is.null( p))
    phi <- parms[ncol( X.p)+1]
  if( !is.null( phi) & is.null( p))
    p <- parms[ncol( X.p)+1]
  
  lambda <- ( mu^( 2-p)) / ( phi*(2-p))
  alpha <- ( 2-p) / ( p-1)
  tau <- phi*(p-1)*mu^(p-1)
  mu.Z <- alpha * tau

  return( -sum( wts * dPoisGam( y, lambda=lambda, mu.Z=mu.Z, alpha=alpha, LOG=TRUE)))
}


"ldTweedie.lp.deriv" <-
function ( parms, y, X.p, offsetty, phi, p, wts=rep( 1, length( y))) 
{
  mu <- exp( X.p %*% parms[1:ncol( X.p)] + offsetty)

  p.flag <- phi.flag <- FALSE
  if( is.null( phi) & is.null( p)){
    p.flag <- phi.flag <- TRUE
    phi <- parms[ncol( X.p) + 1]
    p <- parms[ncol( X.p) + 2]
  }
  if( is.null( phi) & !is.null( p)){
    phi <- parms[ncol( X.p)+1]
    phi.flag <- TRUE
  }
  if( !is.null( phi) & is.null( p)){
    p <- parms[ncol( X.p)+1]
    p.flag <- TRUE
  }

  lambda <- ( mu^( 2-p)) / ( phi*(2-p))
  alpha <- ( 2-p) / ( p-1)
  tau <- phi*(p-1)*mu^(p-1)
  mu.Z <- alpha * tau

  dTweedparms <- -wts * dPoisGamDerivs( y, lambda=lambda, mu.Z=mu.Z, alpha=alpha)

  DTweedparmsDmu <- matrix( c( ( mu^(1-p)) / phi, alpha*phi*( ( p-1)^2)*( mu^(p-2)), rep( 0, length( mu))), nrow=3, byrow=T)
  tmp <- rowSums( dTweedparms * t( DTweedparmsDmu))
  tmp <- tmp * mu
  tmp <- apply( X.p, 2, function( x) x*tmp)
  
  derivs <- colSums( tmp)
  
  if( phi.flag){
    DTweedparmsDphi <- matrix( c( -( ( mu^(2-p)) / ( ( phi^2)*(2-p))), alpha*( p-1)*( mu^( p-1)), rep( 0, length( mu))), nrow=3, byrow=T)
    tmpPhi <- rowSums( dTweedparms * t( DTweedparmsDphi))#vectorised way of doing odd calculation
    derivs <- c( derivs, sum( tmpPhi))
    names( derivs)[length( derivs)] <- "phi"
  }
  if( p.flag){
    dalphadp <- -( 1+alpha) / ( p-1)
    DTweedparmsDp <- matrix( c( lambda*( 1/(2-p) - log( mu)), mu.Z*( dalphadp/alpha + 1/( p-1) + log( mu)), rep( dalphadp, length( y))), nrow=3, byrow=T)
    tmpP <- rowSums( dTweedparms * t( DTweedparmsDp))
    derivs <- c( derivs, sum( tmpP))
    names( derivs)[length( derivs)] <- "p"
  }
  
  return( derivs)
}


"Lmix" <-
function (pars,y,x,G) 
{
   fm <- pars[-1*(1:(G-1))]
    pi <- pars[(1:(G-1))]
    dim(fm) <- c(G,(length(pars)-(G-1))/G)
    pi <- additive.logistic(pi)
    S <- dim(y)[2]
   
   link.fun <- make.link("logit")

    like <- 1
    lf <- matrix(0,G,S)
    lg <- rep(0,S)
    dBi <- array(0,dim=c(S,G,dim(fm)[2]))
    for(s in 1:S){
      tmp.like <- rep(0,G)
      for(g in 1:G){
        lpre <- x%*%fm[g,]
        for(j in 1:dim(fm)[2]){
          dBi[s,g,j] <- sum((y[,s]-link.fun$linkinv(lpre))*x[,j])
        }
        tmp.like[g] <- pi[g]*prod(dbinom(y[,s],1,link.fun$linkinv(lpre),log=F))
        lf[g,s] <- prod(dbinom(y[,s],1,link.fun$linkinv(lpre),log=F))
      }
      lg[s] <- sum(tmp.like)
      like <- like*sum(tmp.like)
    }
    #print(dBi)
  #  print(log(lg))
  #  print(log(lf))
dl.dpi <- rep(0,G)
   for(g in 1:G) dl.dpi[g] <- ( sum( exp(-log(lg) + log(lf[g,]))))
   
    der <- matrix(0,dim(fm)[1],dim(fm)[2])
    for(j in 1:dim(fm)[2]){
      for(g in 1:G){
        for(s in 1:S){
          der[g,j] <- der[g,j]+ exp( -log(lg[s]) + log(pi[g]) + log(lf[g,s])) * dBi[s,g,j]
          #if(g==1 & j == 1) print(c(der[g,j],-log(lg[s]),pi[g],log(lf[g,s]),dBi[s,g,j]))
        }
      }
    }
dpi.deta <- matrix(0,G-1,G)
   ad.trans <- 1+sum(exp(pars[1:(G-1)]))
   for(i in 1:(G-1))
     for(g in 1:G-1){
       if(i==g) {dpi.deta[i,g] <- exp(pars[i])*exp(pars[g])/ad.trans^2}
       else{ dpi.deta[i,g] <- exp(pars[i])/ad.trans - exp(2*pars[g])/ad.trans^2}
     }
   dpi.deta[,G] <- -rowSums(dpi.deta)
   print(dpi.deta)
    print(dl.dpi)
   d1 <- dpi.deta%*%dl.dpi
   print(d1)
    
    list(like,0-der)
}


"logLike.pars" <-
function (pi,coef,sp.form,sp.data,covar.data) 
{

  G <- length(pi)
  pars <- c(additive.logistic(pi,T)[1:(G-1)],unlist(coef))
  
  S <- dim(sp.data)[2]
  if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
  n <- dim(sp.data)[1]
  sp <- rep(sp.name,each=n)
  var.names <- colnames(covar.data)
  covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
  names(covar.data) <- var.names  
  sp.data <- as.data.frame(sp.data)
  data <- data.frame(obs=unlist(sp.data),covar.data)
  names(data)[1] <- as.character(sp.form)[2]
  first.fit <- list(y=data[,1],x=model.matrix(sp.form,data=data))
  logl <- -logLmix(pars,first.fit,G,S,sp,sp.name)
  logl
}



"logLmix" <-
function (pars,first.fit,G,S,sp,spname,out.tau=FALSE) 
{
   tau <- matrix(0,S,G)
##tau,out.tau=FALSE
  if(G>1){
    fm <- pars[-1*(1:(G-1))]
    pi <- pars[(1:(G-1))]
##    fm <- tau[-1*((length(tau)-(G-2)):length(tau))]
    #pi <- tau[((length(tau)-(G-2)):length(tau))]
    dim(fm) <- c(G,(length(pars)-(G-1))/G)
    ##pi[G] <- 1-sum(pi)
    pi <- additive.logistic(pi)
  } else{
    fm <- tau[1:(length(pars)-1)]
    dim(fm) <- c(1,length(fm))
    pi <- 1
  }
  
  link.fun <- make.link("logit")

  
 log.like <- 0
  for(j in 1:S){
    sel.sp <- which(sp==spname[j])
    tmp.like <- rep(0,G)
    for(i in 1:G){
      if(length(fm[i,])==1){lpre <- first.fit$x[sel.sp,]*fm[i,]
      }else{      lpre <- first.fit$x[sel.sp,]%*%fm[i,]}
      tmp.like[i] <- sum(dbinom(first.fit$y[sel.sp],1,link.fun$linkinv(lpre),log=T))
    }
    eps <- max(tmp.like)
    log.like <- log.like +  (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
    tau[j,] <- log(pi) + tmp.like - (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
  }
  
  if(out.tau)return(list(logl=log.like,tau=tau))
  log.like
}


"logLmix.gaussian" <-
function (pars,first.fit,G,S,sp,spname,out.tau=FALSE) 
{
    tau <- matrix(0,S,G)
##tau,out.tau=FALSE
  if(G>1){
    fm <- pars[-1*(1:((G-1)+G))]  ## remove pi
    sp.int <- fm[(length(fm)-(S*G-1)):length(fm)]
    fm <- fm[-1*((length(fm)-(S*G-1)):length(fm))]
    pi <- pars[(1:(G-1))]
    theta <- pars[G:(G+G-1)]
    
##    fm <- tau[-1*((length(tau)-(G-2)):length(tau))]
    #pi <- tau[((length(tau)-(G-2)):length(tau))]
    dim(sp.int) <- c(S,G)
    dim(fm) <- c(G,length(fm)/G)
    ##pi[G] <- 1-sum(pi)
    pi <- additive.logistic(pi)
    
  } else{
    return(0)
    fm <- tau[1:(length(pars)-1)]
    dim(fm) <- c(1,length(fm))
    pi <- 1
  }
  
 
  
 log.like <- 0
  for(j in 1:S){
    sel.sp <- which(sp==spname[j])
    tmp.like <- rep(0,G)
    for(i in 1:G){
      lpre <- cbind(1,first.fit$x[sel.sp,])%*%c(sp.int[j,i],fm[i,])
      ##tmp.like[i] <- sum(dpois(first.fit$y[sel.sp],exp(lpre),log=T))
      tmp.like[i] <- sum(dnorm(first.fit$y[sel.sp],mean=lpre,sd=theta[i],log=T))
      eps <- max(tmp.like)
      log.like <- log.like +  (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
      tau[j,] <- log(pi) + tmp.like - (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
    }
  }
  if(out.tau)return(list(logl=log.like,tau=tau))
  log.like
}


"logLmix.nbinom" <-
function (pars,first.fit,G,S,sp,spname,out.tau=FALSE) 
{
    tau <- matrix(0,S,G)
##tau,out.tau=FALSE
  if(G>1){
    fm <- pars[-1*(1:((G-1)))]  ## remove pi
    ##sp.int <- fm[(length(fm)-(S*G-1)):length(fm)]
    sp.int <- fm[(length(fm)-(2*S-1)):(length(fm)-S)]
    sp.dispersion <- fm[(length(fm)-(S-1)):length(fm)]
    
    ##fm <- fm[-1*((length(fm)-(S*G-1)):length(fm))]
    fm <- fm[-1*((length(fm)-(2*S-1)):length(fm))]
    pi <- pars[(1:(G-1))]
    theta <- pars[G:(G+G-1)]
    
##    fm <- tau[-1*((length(tau)-(G-2)):length(tau))]
    #pi <- tau[((length(tau)-(G-2)):length(tau))]
    ##dim(sp.int) <- c(S,G)
    dim(fm) <- c(G,length(fm)/G)
 
    ##pi[G] <- 1-sum(pi)
    pi <- additive.logistic(pi)
    
  } else{
    return(0)
    fm <- tau[1:(length(pars)-1)]
    dim(fm) <- c(1,length(fm))
    pi <- 1
  }
  
 
  
 log.like <- 0
  for(j in 1:S){
    sel.sp <- which(sp==spname[j])
    tmp.like <- rep(0,G)
    for(i in 1:G){
      lpre <- cbind(1,first.fit$x[sel.sp,])%*%c(sp.int[j],fm[i,])+first.fit$offset[sel.sp]##,i],fm[i,])
      ##tmp.like[i] <- sum(dpois(first.fit$y[sel.sp],exp(lpre),log=T))
      tmp.like[i] <- sum(dnbinom(first.fit$y[sel.sp],mu=exp(lpre),size=sp.dispersion[j],log=T))
    }
   # print(tmp.like)
      eps <- max(tmp.like)
      log.like <- log.like +  (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
      tau[j,] <- log(pi) + tmp.like - (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
    }
  
  if(out.tau)return(list(logl=log.like,tau=tau))
  log.like
}


"logLmix.tweedie" <-
function (pars,first.fit,G,S,sp,spname,out.tau=FALSE) 
{
    tau <- matrix(0,S,G)
##tau,out.tau=FALSE
  if(G>1){
    pi <- pars[which(names(pars)=="pi")]
    phi <- pars[which(names(pars)=="phi")]
    p <- rep(1.6,S)
    fm <- pars[which(names(pars)=="coef")]
    sp.int <- pars[which(names(pars)=="int")]
    ##fm <- pars[-1*(1:((G-1)+G))]  ## remove pi
    ##sp.int <- fm[(length(fm)-(S*G-1)):length(fm)]
    ##fm <- fm[-1*((length(fm)-(S*G-1)):length(fm))]
    ##pi <- pars[(1:(G-1))]
    ##theta <- pars[G:(G+G-1)]
    
##    fm <- tau[-1*((length(tau)-(G-2)):length(tau))]
    #pi <- tau[((length(tau)-(G-2)):length(tau))]
    ##dim(sp.int) <- c(S,G)
    dim(fm) <- c(G,length(fm)/G)
    ##pi[G] <- 1-sum(pi)
    pi <- additive.logistic(pi)
    
  } else{
    return(0)
    fm <- tau[1:(length(pars)-1)]
    dim(fm) <- c(1,length(fm))
    pi <- 1
  }
  
 
  
 log.like <- 0
  for(j in 1:S){
    sel.sp <- which(sp==spname[j])
    tmp.like <- rep(0,G)
    for(i in 1:G){
      lpre <- cbind(1,first.fit$x[sel.sp,])%*%c(sp.int[j],fm[i,])+first.fit$offset[sel.sp]
      ##tmp.like[i] <- sum(dpois(first.fit$y[sel.sp],exp(lpre),log=T))
      tmp.like[i] <- sum(dTweedie(y=first.fit$y[sel.sp],mu=exp(lpre),phi=phi[j],p=p[j],LOG=T))
    }
##print(tmp.like)
    eps <- max(tmp.like)
      log.like <- log.like +  (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
      tau[j,] <- log(pi) + tmp.like - (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
    
  }
  if(out.tau)return(list(logl=log.like,tau=tau))
  log.like
}


"mix.residuals" <-
function (fmM,form,datsp,sp) 
{
   cat("calculating residuals \n")
  link.fun <- make.link("logit")
  x <- model.matrix(form,data=datsp)
  spname <- unique(sp)
  S <- length(spname)
  G <- ncol(fmM$tau)

  PIT <- matrix(NA,S,G)
  for(g in 1:G){
    for(s in 1:length(spname)){
      sel.sp <- which(sp==spname[s])
      t.obs <- sum(datsp$obs[sel.sp])
      pre <- link.fun$linkinv(x[sel.sp,]%*%fmM$coef[g,])
      obs <- datsp$obs[sel.sp]
      dis <- distr.binom(pre)
#      PIT[s,g] <- qnorm(sum(dis[1:(length(dis)-1)],dis[length(dis)]/2))
      nSucc <- sum( obs)
      transfo <- sum( dis[1:nSucc],dis[nSucc+1]/2)
      transfo <- min( transfo, 1)
      ##transfor <- max( transfo, 0)
      PIT[s,g] <- qnorm( transfo)
    }
  }
  PIT 
}


"nd2" <-
function( x0, f, m=NULL, D.accur=4, ...) {
# A function to compute highly accurate first-order derivatives 
# From Fornberg and Sloan (Acta Numerica, 1994, p. 203-267; Table 1, page 213)  
# Adapted by Scott Foster from code nicked off the net 2007

  D.n <- length( x0)
  if ( is.null( m)) {
    D.f0 <- f(x0, ...)
    m <- length( D.f0) 
  }
  if ( D.accur == 2) {
    D.w <- tcrossprod( rep( 1, m),c( -1/2, 1/2))
    D.co <- c( -1, 1) 
  }
  else {
    D.w <- tcrossprod( rep( 1, m),c( 1/12, -2/3, 2/3, -1/12))
    D.co <- c( -2, -1, 1, 2) 
  }
  D.n.c <- length( D.co)
  macheps <- .Machine$double.eps
  D.h <- macheps^( 1/3)*abs( x0)
  D.deriv <- matrix( NA, nrow=m, ncol=D.n)
  for ( ii in 1:D.n) {
    D.temp.f <- matrix( 0, m, D.n.c)
    for ( jj in 1:D.n.c) {
      D.xd <- x0+D.h[ii]*D.co[jj]*( 1:D.n == ii)
      D.temp.f[,jj] <- f( D.xd, ...) 
    }
    D.deriv[,ii] <- rowSums( D.w*D.temp.f)/D.h[ii] 
  }
  return( as.double( D.deriv))
}


"nH2" <-
function( pt, fun, accur=c(4,4), type="H.Diag", ...) {
# A function to compute highly accurate second order Hessian diags and other off-diags
# partially from Fornberg and Sloan (Acta Numerica, 1994, p. 203-267; Table 1, page 213)
# Adapted by Scott Foster from code nicked off the net 2007

  H.n <- length( pt)
  derivs <- function( d.x0, ...) { nd2( x0=d.x0, f=fun, m=1, D.accur=accur[2], ...) }
  Hes <- nd2( x0=pt, f=derivs, D.accur=accur[2], ...)
  Hes <- matrix( Hes, nrow=length( pt))
  Hes <- ( Hes+t( Hes))/2

  if ( type == "H.Diag") {
    macheps <- .Machine$double.eps
    H.h <- macheps^(1/4)*abs(pt)
    H.f0 <- fun( pt, ...)
    H.m <- length( H.f0)
    if ( accur[1] == 2) {
      H.w <- tcrossprod( rep( 1, H.m), c( 1, -2, 1))
      H.co <- c( -1, 0, 1) 
    }
    else {
      H.w <- tcrossprod( rep( 1, H.m), c( -1/12, 4/3, -5/2, 4/3, -1/12))
      H.co <- c( -2, -1, 0, 1, 2) 
    }
    H.n.c <- length( H.co)
    Hes.diag <- double( length=H.n)
    for ( ii in 1:H.n) {
      H.temp.f <- matrix( 0, H.m, H.n.c)
      for ( jj in 1:H.n.c) {
        if ( H.co[jj] != 0) {
          H.xd <- pt+H.h[ii]*H.co[jj]*( 1:H.n == ii)
          H.temp.f[,jj] <- fun( H.xd, ...) 
        }
        else
          H.temp.f[,jj] <- H.f0 
      }
      Hes.diag[ii] <- rowSums( H.w*H.temp.f)/( H.h[ii]^2) 
    } 
    diag( Hes) <- Hes.diag 
  }
  return( Hes)
}


".onLoad" <-
function (libname, pkgname) 
{
  ##require(MASS)
    # Generic DLL loader
    dll.path <- file.path( libname, pkgname, 'libs')
    if( nzchar( subarch <- .Platform$r_arch))
      dll.path <- file.path( dll.path, subarch)
    this.ext <- paste( sub( '.', '[.]', .Platform$dynlib.ext, 
fixed=TRUE), '$', sep='')

    dlls <- dir( dll.path, pattern=this.ext, full.names=FALSE)
    names( dlls) <- dlls
    if( length( dlls))
      lapply( dlls, function( x) library.dynam( sub( this.ext, '', x), 
package=pkgname, lib.loc=libname))
}


"predict.archetype" <-
function (object,new.obs,...)
{
  mixture.model <- object
  if(class(mixture.model)[2]=="bernoulli"){
  G <- length(mixture.model$pi)
  covar <- mixture.model$covar[-(1:(G-1)),-(1:(G-1))]
  coef <- mixture.model$coef
  model.fm <- as.formula(mixture.model$formula)
  model.fm[[2]] <- NULL
  X <- model.matrix(model.fm,new.obs)
  link.fun <- make.link("logit")
  outvar <- matrix(NA,dim(X)[1],G)
  outpred <- matrix(NA,dim(X)[1],G)
  colnames(outvar) <- colnames(outpred) <- paste("G",1:G,sep=".")
  for(g in 1:G){
    lp <- as.numeric(X%*%coef[g,])
    outpred[,g] <- link.fun$linkinv(lp)
    dhdB <- (exp(lp)/(1+exp(lp)))*X - exp(lp)^2/((1+exp(lp))^2)*X
    c2 <- covar[seq(g,dim(covar)[1],G),seq(g,dim(covar)[1],G)]
    for(k in 1:dim(X)[1]){
      outvar[k,g] <- (dhdB[k,]%*%c2)%*%(dhdB[k,])
    }
  }
}
  if(class(mixture.model)[2]=="negbin"|class(mixture.model)[2]=="tweedie"){
    G <- length(mixture.model$pi)
    ##covar <- mixture.model$covar[-1*c(1:(G-1),(dim(mixture.model$covar)[1]-2*length(mixture.model$sp.intercept)+1):dim(mixture.model$covar)[1]),-1*c(1:(G-1),(dim(mixture.model$covar)[1]-2*length(mixture.model$sp.intercept)+1):dim(mixture.model$covar)[1])] # remove pi,sp.int,sp.disp from matrix
    covar <- mixture.model$covar[-1*c(1:(G-1),(dim(mixture.model$covar)[1]-length(mixture.model$sp.intercept)+1):dim(mixture.model$covar)[1]),-1*c(1:(G-1),(dim(mixture.model$covar)[1]-length(mixture.model$sp.intercept)+1):dim(mixture.model$covar)[1])] # remove pi,sp.disp from matrix
    sp.int <- mixture.model$sp.intercept
    coef <- mixture.model$coef
    model.fm <- as.formula(mixture.model$formula)
    model.fm[[2]] <- NULL
    X <- cbind(model.matrix(model.fm,new.obs),1)
    offset <- model.frame(model.fm, data = new.obs)
    offset <- model.offset(offset)
    if(is.null(offset)) offset <-  rep(0,row(X))

    outvar <- matrix(NA,dim(X)[1],G)
    outpred <- matrix(NA,dim(X)[1],G)
    colnames(outvar) <- colnames(outpred) <- paste("G",1:G,sep=".")

    for(g in 1:G){
      s.outvar <- matrix(NA,dim(X)[1],length(sp.int))
      s.outpred <- matrix(NA,dim(X)[1],length(sp.int))
    
      for(s in 1:length(sp.int)){      
        lp <- as.numeric(X%*%c(coef[g,],sp.int[s])+offset)
        s.outpred[,s] <- exp(lp)
        dhdB <- exp(lp)*X
        c2 <- covar[c(seq(g,G*(dim(X)[2]-1),G),G*(dim(X)[2]-1)+s),c(seq(g,G*(dim(X)[2]-1),G),G*(dim(X)[2]-1)+s)]
        for(k in 1:dim(X)[1]){
          s.outvar[k,s] <- (dhdB[k,]%*%c2)%*%(dhdB[k,])
        }
      }
      outpred[,g] <- apply(s.outpred*rep(mixture.model$tau[,g],each=dim(X)[1]),1,mean)/sum(mixture.model$tau[,g])
      outvar[,g]<- apply(s.outvar*rep(mixture.model$tau[,g],each=dim(X)[1]),1,mean)/sum(mixture.model$tau[,g])
    }
  }
  list(fit=outpred,se.fit=sqrt(outvar))
}


"print.archetype" <-
function (x,...) 
{
  cat("\nMixing probabilities\n")
  print(x$pi)
  cat("\nCoefficents\n")
  print(x$coef)
  if(!is.na(x$se[1])){
    cat("\nStandard Errors of coefficents\n")
    print(x$se)
  }
  cat("\nPosterior Probabilities\n")
  print(x$tau)
}


"pTweedie" <-
function ( quant, mu, phi, p) 
{
  lambda <- ( mu^( 2-p)) / ( phi*(2-p))
  alpha <- ( 2-p) / ( p-1)
  tau <- phi*(p-1)*mu^(p-1)
  mu.Z <- alpha * tau
ps <- 0
 # ps <- pPoisGam( quant, lambda, mu.Z, alpha)
 
#  require( tweedie)
#  ps <- ptweedie( as.numeric( q), mu=as.numeric( mu), phi=as.numeric( phi), power=as.numeric( p))
  
  return( ps)
}


"rPoisGam" <-
function ( n, lambda, mu.Z, alpha)
{
  mu.N <- lambda
#simulate n random variables from the same compound poisson distribution
  my.fun <- function (parms)
    return( rgamma( n=1, scale=parms[3], shape=parms[1]*parms[2]))
#    return( sum( rgamma( n=parms[1], scale=parms[3], shape=parms[2])))

  tmp <- c( n, length( mu.N), length(mu.Z), length( alpha))
  names( tmp) <- c( "n","mu.N","mu.Z","alpha")
  if( !all( is.element( tmp[-1], c( 1, tmp[1])))) {
    print( "rPoisGam: error -- length of arguments are not compatible")
    return( tmp)
  }
  if( tmp["mu.N"]==1)
    mu.N <- rep( mu.N, tmp["n"])
  if( tmp["mu.Z"]==1)
    mu.Z <- rep( mu.Z, tmp["n"])
  if( tmp["alpha"]==1)
    alpha <- rep( alpha, tmp["n"])

  np <- matrix( rpois( n, mu.N), ncol=1)
  beta <- mu.Z / alpha
  y <- apply( cbind( np, alpha, beta), 1, my.fun)

  return( y)
}


"rTweedie" <-
function ( n, mu, phi, p)
{
  lambda <- ( mu^( 2-p)) / ( phi*(2-p))
  alpha <- ( 2-p) / ( p-1)
  tau <- phi*(p-1)*mu^(p-1)
  mu.Z <- alpha * tau
 
  rans <- rPoisGam( n, lambda, mu.Z, alpha)
  return( rans)
}


"SpeciesMix" <-
function (sp.form,sp.data,covar.data,G=2, pars=NA, em.prefit=TRUE,em.steps=3, em.refit = 1 , dist="bernoulli",est.var = FALSE,residuals=FALSE,trace=TRUE) 
{
  
  if(dist=="bernoulli") return(SpeciesMix.bernoulli(sp.form,sp.data,covar.data,G, pars, em.prefit,em.steps, em.refit ,est.var,residuals,trace,r1=FALSE))
  ##  if(dist=="bernoulli-random") return(SpeciesMix.bernoulli(sp.form,sp.data,covar.data,G, pars, em.prefit,em.steps, em.refit ,est.var,residuals,trace,r1=TRUE))
  if(dist=="negbin") return(SpeciesMix.nbinom(sp.form,sp.data,covar.data,G, pars, em.prefit,em.steps, em.refit ,est.var,residuals,trace))
  if(dist=="tweedie") return(SpeciesMix.tweedie(sp.form,sp.data,covar.data,G, pars, em.prefit,em.steps, em.refit ,est.var,residuals,trace))
  ##if(dist=="gaussian") return(SpeciesMix.gaussian(sp.form,sp.data,covar.data,G, pars, em.prefit,em.steps, em.refit ,est.var,residuals,trace))
  print("incorrect distribution type, options are bernoulli, negbin or tweedie")
  return(0)
}


"SpeciesMix.bernoulli" <-
function (sp.form,sp.data,covar.data,G=2, pars=NA, em.prefit=TRUE,em.steps=4, em.refit = 1 , est.var = FALSE,residuals=FALSE,trace=TRUE,r1=FALSE) 
{
  t.covar.data <- covar.data
  t.sp.data <- sp.data
  sp.form <- update.formula(sp.form,obs~1+.)
  if(em.prefit | G==1){
    prefit <- SpeciesMix.em(sp.form,sp.data,covar.data,G,ite.max=em.steps,em.refit=em.refit,r1)
    pars <- c(additive.logistic(prefit$pi,T)[1:(G-1)],unlist(prefit$coef))
  }
  S <- dim(sp.data)[2]
  if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
  n <- dim(sp.data)[1]
  sp <- rep(sp.name,each=n)
  if(G==1) return(prefit)
  var.names <- colnames(covar.data)
  covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
  names(covar.data) <- var.names  
  sp.data <- as.data.frame(sp.data)
  data <- data.frame(obs=as.numeric(unlist(sp.data)),covar.data)
  names(data)[1] <- as.character(sp.form)[2]
  fmM.out <- fitmix.cpp(sp.form,data,sp,G,pars=pars,calc.hes=est.var,r1)
  while(fmM.out$logl==0) {
    prefit <- SpeciesMix.em(sp.form,t.sp.data,t.covar.data,G,ite.max=em.steps,em.refit=1,r1)
    pars <- c(additive.logistic(prefit$pi,T)[1:(G-1)],unlist(prefit$coef))
    fmM.out <- fitmix.cpp(sp.form,data,sp,G,pars=pars,calc.hes=est.var,r1)
  }

  rownames(fmM.out$tau) <- sp.name ## add names to taus
  fmM.out$se <- NA
  if(est.var){
    fmM.out$covar <- try(solve(fmM.out$hessian))
     if(class(fmM.out$covar)!="try-error"){
     ##colnames(fmM.out$covar) <- rownames(fmM.out$covar) <- names(fmM.out$gradient)
     tmp <- sqrt(diag(fmM.out$covar))
     tmp <- tmp[(G):length(tmp)]
     fmM.out$se <- matrix(tmp,G,ncol(fmM.out$coef))
     colnames(fmM.out$se) <- colnames(fmM.out$coef)
     rownames(fmM.out$se) <- rownames(fmM.out$coef)
   }
  }
  if(residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
  fmM.out$formula <- sp.form
  class(fmM.out) <- c("archetype","bernoulli")
  fmM.out
}


"SpeciesMix.em" <-
function (sp.form,sp.data,covar.data,G=2,em.refit=1,ite.max=500, est.var = FALSE,residuals=FALSE,trace=TRUE,r1=FALSE) 
{
  S <- dim(sp.data)[2]
  if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
  n <- dim(sp.data)[1]
  sp <- rep(sp.name,each=n)

  var.names <- colnames(covar.data)
  covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
  names(covar.data) <- var.names  
  sp.data <- as.data.frame(sp.data)
  data <- data.frame(obs=unlist(sp.data),covar.data)
  
  fmM.out <- fitMix(sp.form,data,sp,G,ite.max,trace,r1)
  if(em.refit>1)
    for(i in 2:em.refit){
      fmM <- fitMix(sp.form,data,sp,G,ite.max,trace,r1)
      if(fmM$logl>fmM.out$logl) fmM.out <- fmM
    }
  
  if(est.var){
    var <- 0
    t.pi <- additive.logistic(fmM.out$pi,TRUE)
    parms <- c(t.pi[1:(G-1)],unlist(fmM.out$coef))
    first.fit <- list(y=data[,1],x=model.matrix(sp.form,data=data))
    fun.est.var <- function(x){-logLmix(x,first.fit,G,S,sp,sp.name)}
    var <- solve( nH2( pt=parms, fun=fun.est.var))
    colnames( var) <- rownames( var) <- names( parms)
    fmM.out$covar <- var
  }
  if(residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
  fmM.out$formula <- sp.form

  fmM.out
}


"SpeciesMix.em.gaussian" <-
function (sp.form,sp.data,covar.data,G=2,em.refit=1,ite.max=500, est.var = FALSE,residuals=FALSE,trace=TRUE) 
{
  S <- dim(sp.data)[2]
  if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
  n <- dim(sp.data)[1]
  sp <- rep(sp.name,each=n)

  var.names <- colnames(covar.data)
  covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
  names(covar.data) <- var.names  
  sp.data <- as.data.frame(sp.data)
  data <- data.frame(obs=unlist(sp.data),covar.data)
  
  fmM.out <- fitmix.gaussian(sp.form,data,sp,G,ite.max,trace)
  if(em.refit>1)
    for(i in 2:em.refit){
      fmM <- fitmix.gaussian(sp.form,data,sp,G,ite.max,trace)
      if(fmM$logl>fmM.out$logl) fmM.out <- fmM
    }
  ##est.var <- FALSE
  if(est.var){
    var <- 0
    t.pi <- additive.logistic(fmM.out$pi,TRUE)
    ##parms <- c(t.pi[1:(G-1)],unlist(fmM.out$coef))
    parms <- c(t.pi[1:(G-1)],fmM.out$theta,unlist(fmM.out$coef),unlist(fmM.out$sp.intercept))
    first.fit <- list(y=data[,1],x=model.matrix(sp.form,data=data)[,-1])
    fun.est.var <- function(x){-logLmix.gaussian(x,first.fit,G,S,sp,sp.name)}
    var <- solve( nH2( pt=parms, fun=fun.est.var))
    colnames( var) <- rownames( var) <- names( parms)
    fmM.out$covar <- var
  }
  residuals <- FALSE
  if(residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
  fmM.out$formula <- sp.form

  fmM.out
}


"SpeciesMix.em.nbinom" <-
function (sp.form,sp.data,covar.data,G=2,em.refit=1,ite.max=500, est.var = FALSE,residuals=FALSE,trace=TRUE) 
{
  S <- dim(sp.data)[2]
  if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
  n <- dim(sp.data)[1]
  sp <- rep(sp.name,each=n)

  var.names <- colnames(covar.data)
  covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
  names(covar.data) <- var.names  
  sp.data <- as.data.frame(sp.data)
  data <- data.frame(obs=unlist(sp.data),covar.data)
  
  fmM.out <- fitmix.nbinom(sp.form,data,sp,G,ite.max,trace)
  if(em.refit>1)
    for(i in 2:em.refit){
      fmM <- fitmix.nbinom(sp.form,data,sp,G,ite.max,trace)
      if(fmM$logl>fmM.out$logl) fmM.out <- fmM
    }
  ##est.var <- FALSE
  if(est.var){
    var <- 0
    t.pi <- additive.logistic(fmM.out$pi,TRUE)
    ##parms <- c(t.pi[1:(G-1)],unlist(fmM.out$coef))
    parms <- c(t.pi[1:(G-1)],unlist(fmM.out$coef),fmM.out$sp.intercept,rep(1,S))##unlist(fmM.out$sp.intercept))
     first.fit <- list(y=data[,1],x=model.matrix(sp.form,data=data)[,-1])
   # first.fit$x <- cbind(sp.mat,first.fit$x)
    fun.est.var <- function(x){-logLmix.nbinom(x,first.fit,G,S,sp,sp.name)}
    #
    deriv <- nd2(parms,fun.est.var)
    var <- solve( nH2( pt=parms, fun=fun.est.var))
   colnames( var) <- rownames( var) <- names( parms)
   fmM.out$covar <- var
  }
  residuals <- FALSE
  if(residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
  fmM.out$formula <- sp.form

  fmM.out
}


"SpeciesMix.em.tweedie" <-
function (sp.form,sp.data,covar.data,G=2,em.refit=1,ite.max=500, est.var = FALSE,residuals=FALSE,trace=TRUE) 
{
  S <- dim(sp.data)[2]
  if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
  n <- dim(sp.data)[1]
  sp <- rep(sp.name,each=n)

  var.names <- colnames(covar.data)
  covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
  names(covar.data) <- var.names  
  sp.data <- as.data.frame(sp.data)
  data <- data.frame(obs=unlist(sp.data),covar.data)
  
  fmM.out <- fitmix.tweedie(sp.form,data,sp,G,ite.max,trace)
  if(em.refit>1)
    for(i in 2:em.refit){
      fmM <- fitmix.tweedie(sp.form,data,sp,G,ite.max,trace)
      if(fmM$logl>fmM.out$logl) fmM.out <- fmM
    }
  est.var <- FALSE
  if(est.var){
    var <- 0
    t.pi <- additive.logistic(fmM.out$pi,TRUE)
    ##parms <- c(t.pi[1:(G-1)],unlist(fmM.out$coef))
    parms <- c(t.pi[1:(G-1)],unlist(fmM.out$coef),unlist(fmM.out$sp.intercept),fmM.out$theta)
    first.fit <- list(y=data[,1],x=model.matrix(sp.form,data=data)[,-1])
    fun.est.var <- function(x){-logLmix.tweedie(x,first.fit,G,S,sp,sp.name)}
    var <- solve( nH2( pt=parms, fun=fun.est.var))
    colnames( var) <- rownames( var) <- names( parms)
    fmM.out$covar <- var
  }
  residuals <- FALSE
  if(residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
  fmM.out$formula <- sp.form

  fmM.out
}


"SpeciesMix.gaussian" <-
function (sp.form,sp.data,covar.data,G=2, pars=NA, em.prefit=TRUE,em.steps=4, em.refit = 1 , est.var = FALSE,residuals=FALSE,trace=TRUE) 
{
  t.covar.data <- covar.data
  t.sp.data <- sp.data
  sp.form <- update.formula(sp.form,obs~1+.)
  em.steps=100
  ##if(em.prefit | G==1){
    prefit <- SpeciesMix.em.gaussian(sp.form,sp.data,covar.data,G,ite.max=em.steps,est.var=est.var,em.refit=em.refit)
  return(prefit)
    ##pars <- c(additive.logistic(prefit$pi,T)[1:(G-1)],unlist(prefit$coef))
  ##}
  ##S <- dim(sp.data)[2]
  ##if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
  ##n <- dim(sp.data)[1]
  ##sp <- rep(sp.name,each=n)
  ##if(G==1) return(prefit)
  ##var.names <- colnames(covar.data)
  ##covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
  ##names(covar.data) <- var.names  
  ##sp.data <- as.data.frame(sp.data)
  ##data <- data.frame(obs=as.numeric(unlist(sp.data)),covar.data)
  ##names(data)[1] <- as.character(sp.form)[2]
  ##fmM.out <- fitmix.cpp(sp.form,data,sp,G,pars=pars,calc.hes=est.var)
  ##while(fmM.out$logl==0) {
   ## prefit <- SpeciesMix.em(sp.form,t.sp.data,t.covar.data,G,ite.max=em.steps,em.refit=1)
   ## pars <- c(additive.logistic(prefit$pi,T)[1:(G-1)],unlist(prefit$coef))
    ##fmM.out <- fitmix.cpp(sp.form,data,sp,G,pars=pars,calc.hes=est.var)
  ##}

  ##rownames(fmM.out$tau) <- sp.name ## add names to taus

  ##if(est.var){
  ## fmM.out$covar <- try(solve(fmM.out$hessian))
  ##}
 ## if(residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
  ##fmM.out$formula <- sp.form
  ##class(fmM.out) <- "archetype"
  ##fmM.out
}


"SpeciesMix.nbinom" <-
function (sp.form,sp.data,covar.data,G=2, pars=NA, em.prefit=TRUE,em.steps=3, em.refit = 1 , est.var = FALSE,residuals=FALSE,trace=TRUE) 
{
  t.covar.data <- covar.data
  t.sp.data <- sp.data
  sp.form <- update.formula(sp.form,obs~1+.)
  #em.steps=100
  if(em.prefit | G==1){
    prefit <- SpeciesMix.em.nbinom(sp.form,sp.data,covar.data,G,ite.max=em.steps,est.var=FALSE,em.refit=em.refit)
    ##return(prefit)
    pars <- prefit$pars#c(additive.logistic(prefit$pi,T)[1:(G-1)],unlist(prefit$coef))
  }
  S <- dim(sp.data)[2]
  if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
  n <- dim(sp.data)[1]
  sp <- rep(sp.name,each=n)
  if(G==1) return(prefit)
  var.names <- colnames(covar.data)
  covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
  names(covar.data) <- var.names  
  sp.data <- as.data.frame(sp.data)
  data <- data.frame(obs=as.numeric(unlist(sp.data)),covar.data)
  names(data)[1] <- as.character(sp.form)[2]
  sp.form <- update.formula(sp.form,obs~-1+.)
  fmM.out <- fitmix.nbinom.cpp(sp.form,data,sp,G,pars=pars,calc.hes=est.var)
  #while(fmM.out$logl==0) {
  #  prefit <- SpeciesMix.em.nbinom(sp.form,sp.data,covar.data,G,ite.max=em.steps,est.var=est.var,em.refit=em.refit)
  # pars <- prefit$pars
  # fmM.out <- fitmix.nbinom.cpp(sp.form,data,sp,G,pars=pars,calc.hes=est.var)
  #}

  rownames(fmM.out$tau) <- sp.name ## add names to taus
  fmM.out$se <- NA
  if(est.var){
      fmM.out$covar <- try(solve(fmM.out$hessian))
   if(class(fmM.out$covar)!="try-error"){
     #colnames(fmM.out$covar) <- rownames(fmM.out$covar) <- names(fmM.out$gradient)
     tmp <- sqrt(diag(fmM.out$covar))
     tmp <- tmp[-1*(1:((G-1)))]
     tmp <- tmp[-1*((length(tmp)-(2*S-1)):length(tmp))]
     fmM.out$se <- matrix(tmp,G,ncol(fmM.out$coef))
     colnames(fmM.out$se) <- colnames(fmM.out$coef)
     rownames(fmM.out$se) <- rownames(fmM.out$coef)
   }
      ##t1 <- sqrt(diag(fmM.out$covar))
    ##t1 <- t1[-(1:G-1)]
    ##fmM.out$se.coef <- t1[1:length(fmM.out$coef)]
    ##t1 <- t1[-(1:length(fmM.out$coef))]
    ##dim(fmM.out$se.coef) <- dim(fmM.out$coef)
    
    ##fmM.out$se.int <- t1[1:S]
    ##fmM.out$se.disp <- t1[-(1:S)]
  }
  if(residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
  fmM.out$formula <- sp.form
  class(fmM.out) <- c("archetype","negbin")
  fmM.out
}


"SpeciesMix.tweedie" <-
function (sp.form,sp.data,covar.data,G=2, pars=NA, em.prefit=TRUE,em.steps=3, em.refit = 1 , est.var = FALSE,residuals=FALSE,trace=TRUE) 
{
  t.covar.data <- covar.data
  t.sp.data <- sp.data
  sp.form <- update.formula(sp.form,obs~1+.)
  pars <- NA
  if(em.prefit | G==1){
    prefit <- SpeciesMix.em.tweedie(sp.form,sp.data,covar.data,G,ite.max=em.steps,est.var=est.var,em.refit=em.refit)
  ##return(prefit)
    pars <- prefit$pars##pars <- c(additive.logistic(prefit$pi,T)[1:(G-1)],unlist(prefit$coef))
  }
  S <- dim(sp.data)[2]
  if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
  n <- dim(sp.data)[1]
  sp <- rep(sp.name,each=n)
  if(G==1) return(prefit)
  var.names <- colnames(covar.data)
  covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
  names(covar.data) <- var.names  
  sp.data <- as.data.frame(sp.data)
  data <- data.frame(obs=as.numeric(unlist(sp.data)),covar.data)
  names(data)[1] <- as.character(sp.form)[2]
  sp.form <- update.formula(sp.form,obs~-1+.)
  fmM.out <- fitmix.tweedie.cpp(sp.form,data,sp,G,pars=pars,calc.hes=est.var)
  ##while(fmM.out$logl==0) {
   ## prefit <- SpeciesMix.em(sp.form,t.sp.data,t.covar.data,G,ite.max=em.steps,em.refit=1)
   ## pars <- c(additive.logistic(prefit$pi,T)[1:(G-1)],unlist(prefit$coef))
    ##fmM.out <- fitmix.cpp(sp.form,data,sp,G,pars=pars,calc.hes=est.var)
  ##}

  rownames(fmM.out$tau) <- sp.name ## add names to taus
  fmM.out$se <- NA
  if(est.var){
   ##fmM.out$covar <- try(solve(fmM.out$hessian))
   fmM.out$covar <- try(solve(fmM.out$hessian))
   if(class(fmM.out$covar)!="try-error"){
     colnames(fmM.out$covar) <- rownames(fmM.out$covar) <- names(fmM.out$gradient)
     tmp <- sqrt(diag(fmM.out$covar))
     tmp <- tmp[-1*(1:((G-1)))]
     tmp <- tmp[-1*((length(tmp)-(2*S-1)):length(tmp))]
     fmM.out$se <- matrix(tmp,G,ncol(fmM.out$coef))
     colnames(fmM.out$se) <- colnames(fmM.out$coef)
     rownames(fmM.out$se) <- rownames(fmM.out$coef)
   }
  }
  if(residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
 fmM.out$formula <- sp.form
 class(fmM.out) <- c("archetype","tweedie")
 fmM.out
}


"tglm" <-
function ( mean.form, data, wts=NULL, phi=NULL, p=NULL, inits=NULL, vcov=TRUE, residuals=TRUE, trace=1, iter.max=150) 
{
  if( is.null( wts))
    wts <- rep( 1, nrow( data))

  e<-new.env()
  e$wts <- wts
  environment(mean.form) <- e
  temp.p <- model.frame( mean.form, data=as.data.frame( data), weights=wts)
  y <- model.response( temp.p)
  names( y) <- NULL
  X.p <- model.matrix( mean.form, data)
  offset.p <- model.offset( temp.p)
  if ( is.null( offset.p))
    offset.p <- rep( 0, nrow( temp.p))
  wts1 <- as.vector( model.weights( temp.p))

  fm1 <- NULL
  if( is.null( inits)){
    inits <- rep( 0, ncol( X.p))
    if( trace!=0)
      print( "Obtaining initial mean values from log-linear Poisson model -- this might be stupid")
    abit <- 1
    ystar <- round( y, digits=0)
    fm1 <- glm( ystar~-1+X.p, family=poisson( link="log"), weights=wts1)
    inits <- fm1$coef
  }

  if( is.null( phi) & length( inits)==ncol( X.p)){
    if( is.null( fm1))
      inits <- c( inits, 1)
    else{
      if( trace!=0)
        print( "Obtaining initial dispersion from the smaller of the Pearson or Deviance estimator (or 25)")
      if( is.null( p) & length( inits)==ncol( X.p))
        ptemp <- 1.9
      else
        ptemp <- p
      disDev <- fm1$deviance/( length( y) - length( fm1$coef))
      disPear <- sum( ( wts1 * ( y-fm1$fitted)^2) / ( fm1$fitted^ptemp)) / ( length( y) - length( fm1$coef))
      dis <- min( disDev, disPear, 25)
      inits <- c( inits, dis)
    }
  }

  if( is.null( p) & length( inits)==ncol( X.p) + is.null( phi))
    inits <- c( inits, 1.6)

  if( length( inits) != ncol( X.p) + is.null( phi) + is.null( p)) {
    print( "Initial values supplied are of the wrong length -- please check")
    tmp <- c( length( inits), ncol( X.p) + is.null( phi) + is.null( p))
    names( tmp) <- c( "inits", "nParams")
    return( tmp)
  }

  fmTGLM <- tglm.fit( x=X.p, y=y, wts=wts1, offset=offset.p, inits=inits, phi=phi, p=p, vcov=vcov, residuals=residuals, trace=trace, iter.max=iter.max)

  return( fmTGLM)
}


"tglm.fit" <-
function ( x, y, wts=NULL, offset=rep( 0, length( y)), inits, phi=NULL, p=NULL, vcov=TRUE, residuals=TRUE, trace=1, iter.max=150)
{
  if( trace!=0){
    print( "Estimating parameters")
    if( is.null( phi) & is.null( p))
      cat("iter:", "-logl", colnames(x), "phi", "p", "\n", sep = "\t") 
    if( is.null( phi) & !is.null( p))
      cat("iter:", "-logl", colnames(x), "phi", "\n", sep = "\t") 
    if( !is.null( phi) & is.null( p))
      cat("iter:", "-logl", colnames(x), "p", "\n", sep = "\t") 
    if( !is.null( phi) & !is.null( p))
      cat("iter:", "-logl", colnames(x), "\n", sep = "\t") 
  }

  eps <- 1e-5
  my.lower <- rep( -Inf, ncol( x))
  my.upper <- rep( Inf, ncol( x))
  if( is.null( phi))
    {my.lower <- c( my.lower, eps); my.upper <- c( my.upper, Inf)}
  if( is.null( p))
    {my.lower <- c( my.lower, 1+eps); my.upper <- c( my.upper, 2-eps)}
  
  fm <- nlminb( start=inits, objective=ldTweedie.lp, gradient=ldTweedie.lp.deriv, hessian=NULL, lower=my.lower, upper=my.upper, y=y, X.p=x, offsetty=offset, phi=phi, p=p, control=list(trace=trace, iter.max=iter.max), wts=wts)
    
  parms <- fm$par
  tmp <- colnames( x)
  if( !is.null( phi) & !is.null( p))
    tmp <- tmp
  if( !is.null( phi) & is.null( p))
    tmp <- c( tmp, "p")
  if( is.null( phi) & !is.null( p))
    tmp <- c( tmp, "phi")
  if( is.null( phi) & is.null( p))
    tmp <- c( tmp, "phi", "p")
   
  names( parms) <- tmp
  
  if( vcov){
    if( trace!=0)
      print( "Calculating variance matrix of estimates")
    vcovar <- nd2(x0=parms, f=ldTweedie.lp.deriv, y=y, X.p=x, offsetty=offset, phi=phi, p=p, wts=wts)
    vcovar <- 0.5 * ( vcovar + t( vcovar))
    vcovar <- solve( vcovar)
    rownames( vcovar) <- colnames( vcovar) <- names( parms)
  }
  else{
    if( trace!=0)
      print( "Not calculating variance matrix of estimates")
    vcovar <- NULL
  }

  scores <- -ldTweedie.lp.deriv( parms=parms, y=y, X.p=x, offsetty=offset, phi=phi, p=p, wts=wts) 

  if( trace !=0)
    print( "Calculating means")
  mu <- exp( x %*% parms[1:ncol( x)] + offset)

  if( residuals){
    if( trace!=0)
      print( "Calculating quantile residuals")
    if( is.null( phi)){
      phi1 <- parms[ncol( x)+1]
      if( is.null( p))
        p1 <- parms[ncol(x)+2]
      else
        p1 <- p
    }
    else{
      phi1 <- phi
      if( is.null( p))
        p1 <- parms[ncol( x)+1]
      else
        p1 <- p
    }
    resids <- matrix( rep( pTweedie( y, mu, phi1, p1), 2), ncol=2)
    resids[y==0,1] <- 0.5 * dTweedie( y[y==0], mu[y==0], phi1, p1, LOG=FALSE)
    nzero <- sum( y==0)
    resids[y==0,2] <- runif( nzero, min=rep( 0, nzero), max=2*resids[y==0,1])
    resids <- qnorm( resids)
    colnames( resids) <- c("expect","random")
  }
  else{
    if( trace!=0)
      print( "Not calculating quantile residuals")
    resids <- NULL
  }

  if( trace!=0)
    print( "Done")

  ICadj <- 0
  if( !is.null( phi))
    ICadj <- ICadj+1
  if( !is.null( p))
    ICadj <- ICadj+1 

  AIC <- -2*(-fm$objective) + 2*(length( parms)-ICadj)
  BIC <- -2*(-fm$objective) + log( nrow( x))*(length( parms)-ICadj)

  res <- list( coef=parms, logl=-fm$objective, scores=scores, vcov=vcovar, conv=fm$convergence, message=fm$message, niter=fm$iterations, evals=fm$evaluations, call=match.call(), fitted=mu, residuals=resids, AIC=AIC, BIC=BIC)

  class( res) <- "tglm"

  return( res)


}


"weighted.glm" <-
function (g,first.fit,tau,n,fmM,sp) 
{
  dat.tau <- rep(tau[,g],each=n)
  ##lpre <- first.fit$x%*%fmM[[g]]$coef
  ##f.mix <- glm.fit(x=first.fit$x,y=first.fit$y,weights=dat.tau[,g],family=binomial(),etastart=fmM[[g]]$linear.predictors)
  
  f.mix <- glm.fit(x=first.fit$x,y=first.fit$y,weights=dat.tau,family=binomial(),start=fmM[[g]]$coef)
  return(list(coef=f.mix$coef))  
 
  
  ##list(residuals=f.mix$residuals,fitted=f.mix$fitted,linear.predictors=f.mix$linear.predictors,coef=f.mix$coef)
  ##list(linear.predictors=f.mix$linear.predictors,coef=f.mix$coef)
 
  
}


"weighted.glm.gaussian" <-
function (g,first.fit,tau,n,fmM,sp) 
{
  dat.tau <- rep(tau[,g],each=n)
  sp.name <- unique(sp)
  X <- first.fit$x
  sp.mat <- matrix(0,dim(X)[1],length(sp.name))
    
  for(i in 1:length(sp.name)){
    sp.mat[sp==sp.name[i],i] <- 1
  }
  X <- cbind(sp.mat,X)
  f.mix <- glm.fit(x=X,y=first.fit$y,weights=dat.tau,family=gaussian())
  #f.mix$theta <- 1
 ## f.mix <- glm.fit.nbinom(x=X,y=first.fit$y,weights=dat.tau)
  sp.intercept <- f.mix$coef[1:length(sp.name)]
  sp.intercept[is.na(sp.intercept)] <- 0
  return(list(coef=f.mix$coef[-1:-length(sp.name)],theta=sqrt(f.mix$deviance/f.mix$df.residual),sp.intercept=sp.intercept,fitted=f.mix$fitted))#,lpre=f.mix$linear.predictors))
  
}


"weighted.glm.nbinom" <-
function (g,first.fit,tau,n,fmM,sp) 
{
  dat.tau <- rep(tau[,g],each=n)
  sp.name <- unique(sp)
  X <- first.fit$x
  sp.mat <- matrix(0,dim(X)[1],length(sp.name))
    
  for(i in 1:length(sp.name)){
    sp.mat[sp==sp.name[i],i] <- 1
  }
  X <- cbind(sp.mat,X)
  #f.mix <- glm.fit(x=X,y=first.fit$y,weights=dat.tau,family=poisson())
  #f.mix$theta <- 1
  f.mix <- glm.fit.nbinom(x=X,y=first.fit$y,offset=first.fit$offset,weights=dat.tau)
  sp.intercept <- f.mix$coef[1:length(sp.name)]
  sp.intercept[is.na(sp.intercept)] <- 0
  return(list(coef=f.mix$coef[-1:-length(sp.name)],theta=f.mix$theta,sp.intercept=sp.intercept,fitted=f.mix$fitted))#,lpre=f.mix$linear.predictors))
  
}


"weighted.glm.tweedie" <-
function (g, first.fit, tau, n, fmM, sp) 
{
    dat.tau <- rep(tau[, g], each = n)
    sp.int.offset <- rep(fmM[[g]]$sp.intercept, each = n)+first.fit$offset
    sp.name <- unique(sp)
    X <- first.fit$x
    f.mix <- tglm.fit(x = X, y = first.fit$y, wts = dat.tau, 
        offset = sp.int.offset, vcov = FALSE, residuals = FALSE, 
        trace = 0, inits = fmM[[g]]$coef, phi = fmM[[g]]$phi, 
        p = 1.6)
    return(list(coef = f.mix$coef, phi = fmM[[g]]$phi, p = 1.6, 
        sp.intercept = fmM[[g]]$sp.intercept, fitted = f.mix$fitted))
}

# MVB's workaround for futile CRAN 'no visible blah' check:
globalVariables( package="SpeciesMix",
  names=c( ".Traceback"
    ,"inv"
    ,"x"
    ,"x.t"
    ,"dat.tau"
    ,"tau"
    ,"i"
    ,"n"
    ,"form"
    ,"datsp"
    ,"y"
    ,"obs"
    ,"f.mix"
    ,"sp.int"
    ,"theta"
    ,"X"
    ,"data"
    ,"out"
    ,"S"
    ,"k"
    ,"group"
    ,"s"
    ,"g"
    ,"lgtp"
    ,"p"
    ,"tmp"
    ,"my.fun"
    ,"sp.data"
    ,"covar.data"
    ,"em.prefit"
    ,"em.steps"
    ,"em.refit"
    ,"est.var"
    ,"G"
    ,"sp.form"
    ,"aic"
    ,"bic"
    ,"fm"
    ,"logl"
    ,"covar"
    ,"fmM"
    ,"first.fit"
    ,"MM"
    ,"tol"
    ,"sp.name"
    ,"sp"
    ,"starting.fitem"
    ,"all.betas"
    ,"j"
    ,"fit"
    ,"sp.intercepts"
    ,"fmmvnorm"
    ,"centers"
    ,"B"
    ,"cluster"
    ,"est.tau"
    ,"max.newTau"
    ,"alpha"
    ,"newTau"
    ,"phi"
    ,"new.dist"
    ,"old.dist"
    ,"ii"
    ,"jj"
    ,"mu.N"
    ,"lambda"
    ,"mu.Z"
    ,"res"
    ,"LOG"
    ,"do.checks"
    ,"mu"
    ,"dens"
    ,"tmp.like"
    ,"sel.sp"
    ,"spname"
    ,"link.fun"
    ,"lpre"
    ,"linkinv"
    ,"eps"
    ,"sum.like"
    ,"temp.warn"
    ,"ite"
    ,"logL"
    ,"old.logL"
    ,"t1"
    ,"ite.max"
    ,"fm.out"
    ,"EN"
    ,"d"
    ,"t.pi"
    ,"parms"
    ,"logL.full"
    ,"full.model"
    ,"pars"
    ,"gradient"
    ,"model.type"
    ,"loglike"
    ,"calc.deriv"
    ,"ll"
    ,"hes"
    ,"calc.hes"
    ,"r.logl"
    ,"int.out"
    ,"fm.theta"
    ,"sp.intercept"
    ,"sp.dispersion"
    ,"r.deriv"
    ,"fm.phi"
    ,"fm.p"
    ,"sp.phi"
    ,"sp.p"
    ,"pars.og"
    ,"mustart"
    ,"se"
    ,"X.p"
    ,"offsetty"
    ,"wts"
    ,"p.flag"
    ,"phi.flag"
    ,"dTweedparms"
    ,"DTweedparmsDmu"
    ,"derivs"
    ,"DTweedparmsDphi"
    ,"tmpPhi"
    ,"dalphadp"
    ,"DTweedparmsDp"
    ,"tmpP"
    ,"like"
    ,"lf"
    ,"lg"
    ,"dBi"
    ,"dl.dpi"
    ,"der"
    ,"dpi.deta"
    ,"ad.trans"
    ,"d1"
    ,"var.names"
    ,"log.like"
    ,"out.tau"
    ,"PIT"
    ,"t.obs"
    ,"pre"
    ,"dis"
    ,"nSucc"
    ,"transfo"
    ,"D.n"
    ,"x0"
    ,"m"
    ,"D.f0"
    ,"f"
    ,"..."
    ,"D.accur"
    ,"D.w"
    ,"D.co"
    ,"D.n.c"
    ,"macheps"
    ,"double.eps"
    ,"D.h"
    ,"D.deriv"
    ,"D.temp.f"
    ,"D.xd"
    ,"H.n"
    ,"d.x0"
    ,"fun"
    ,"accur"
    ,"Hes"
    ,"type"
    ,"H.h"
    ,"H.f0"
    ,"H.m"
    ,"H.w"
    ,"H.co"
    ,"H.n.c"
    ,"Hes.diag"
    ,"H.temp.f"
    ,"H.xd"
    ,"dll.path"
    ,"libname"
    ,"pkgname"
    ,"subarch"
    ,"r_arch"
    ,"this.ext"
    ,"dynlib.ext"
    ,"dlls"
    ,"mixture.model"
    ,"object"
    ,"model.fm"
    ,"new.obs"
    ,"outvar"
    ,"outpred"
    ,"lp"
    ,"dhdB"
    ,"c2"
    ,"s.outvar"
    ,"s.outpred"
    ,"ps"
    ,"np"
    ,"rans"
    ,"t.covar.data"
    ,"t.sp.data"
    ,"."
    ,"prefit"
    ,"r1"
    ,"fmM.out"
    ,"fun.est.var"
    ,"e"
    ,"mean.form"
    ,"temp.p"
    ,"offset.p"
    ,"wts1"
    ,"fm1"
    ,"inits"
    ,"abit"
    ,"ystar"
    ,"ptemp"
    ,"disDev"
    ,"disPear"
    ,"fmTGLM"
    ,"iter.max"
    ,"my.lower"
    ,"my.upper"
    ,"par"
    ,"vcovar"
    ,"scores"
    ,"phi1"
    ,"p1"
    ,"resids"
    ,"nzero"
    ,"ICadj"
    ,"objective"
    ,"convergence"
    ,"iterations"
    ,"evaluations"
    ,"sp.mat"
    ,"sp.int.offset"
))

