# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Copyright - Maurice Berk (maurice@mauriceberk.com) (http://www2.imperial.ac.uk/~mab201)
#
# block diagonal matrix code thanks to http://tolstoy.newcastle.edu.au/R/help/04/05/1320.html

sme <- function(object,tme,ind,verbose=F,lambda.mu=NULL,lambda.v=NULL,maxIter=500,knots=NULL,zeroIntercept=F,deltaEM=1e-3,deltaNM=1e-3,criteria="AICc",...)
{
  UseMethod("sme")
}

sme.default <- function(object,tme,ind,verbose=F,lambda.mu=NULL,lambda.v=NULL,maxIter=500,knots=NULL,zeroIntercept=F,deltaEM=1e-3,deltaNM=1e-3,criteria="AICc",...)
{
  y <- object
  ind <- as.factor(ind)
  ind.unique <- sort(unique(ind))
  y <- y[order(ind)]
  tme <- tme[order(ind)]
  ind <- ind[order(ind)]

  call <- match.call()

  if(!is.null(knots))
  {
    require(splines)
  
    max.tme <- max(tme)
    min.tme <- min(tme)
    B <- ns(c(min.tme,knots,max.tme),knots=knots,intercept=T)

    G <- roughnessMatrix(incidenceMatrix(c(min.tme,knots,max.tme)))
    G <- t(B) %*% G %*% B
  }
  else
  {
    G <- roughnessMatrix(incidenceMatrix(tme))
  }

  yi <- split(y,ind)
  tmei <- split(tme,ind)
  Ni <- sapply(yi,length)

  if(is.null(knots))
  {
    X <- incidenceMatrix(tme)
    Xi <- split.data.frame(X,ind)
  }
  else
  {
    X <- ns(tme,knots=knots,intercept=T)
    Xi <- split.data.frame(X,ind)
  }

  if(is.null(lambda.mu) || is.null(lambda.v))
  {
    res.em <- .C("SMEOptimization",
                 y=as.double(y),
                 tme=as.double(tme),
                 ind=as.integer(ind),
                 X=as.double(X),
                 N=as.integer(length(y)),
                 n=as.integer(length(yi)),
                 Ni=as.integer(sapply(yi,length)),
                 p=as.integer(ncol(X)),
                 lambdaMu=double(1),
                 lambdaV=double(1),
                 G=as.double(G),
                 mu=double(ncol(X)),
                 sigmaSquared=double(1),
                 D=double(ncol(X) * ncol(X)),
                 v=double(ncol(X) * length(yi)),
                 likelihood=double(1),
                 dfMu=double(1),
                 dfV=double(1),
                 zeroIntercept=as.integer(zeroIntercept),
                 iterations=integer(1),
                 maxIterations=as.integer(maxIter),
                 deltaEM=as.double(deltaEM),
                 deltaNM=as.double(deltaNM),
                 criteria=as.integer(match(criteria,c("AIC","AICc","BICN","BICn"))),
                 verbose=as.integer(verbose),
                 info=integer(1),
                 PACKAGE="sme")
  }
  else
  {
    res.em <- .C("SME",
                 y=as.double(y),
                 tme=as.double(tme),
                 ind=as.integer(ind),
                 X=as.double(X),
                 N=as.integer(length(y)),
                 n=as.integer(length(yi)),
                 Ni=as.integer(sapply(yi,length)),
                 p=as.integer(ncol(X)),
                 lambdaMu=as.double(lambda.mu),
                 lambdaV=as.double(lambda.v),
                 G=as.double(G),
                 mu=double(ncol(X)),
                 sigmaSquared=double(1),
                 D=double(ncol(X) * ncol(X)),
                 v=double(ncol(X) * length(yi)),
                 likelihood=double(1),
                 dfMu=double(1),
                 dfV=double(1),
                 zeroIntercept=as.integer(zeroIntercept),
                 iterations=integer(1),
                 maxIterations=as.integer(maxIter),
                 deltaEM=as.double(deltaEM),
                 verbose=as.integer(verbose),
                 info=integer(1),
                 PACKAGE="sme")
  }

  return.value <- list()
  if(is.null(knots))
  {
    return.value$coefficients <- rbind(as.vector(res.em$mu),t(matrix(res.em$v,ncol=res.em$n)))
    colnames(return.value$coefficients) <- attr(X,"tau")
  }
  else
  {
    return.value$coefficients <- rbind(as.vector(B %*% res.em$mu),do.call(rbind,lapply(split(res.em$v,rep(1:res.em$n,each=res.em$p)),function(v) as.vector(B %*% v))))
    colnames(return.value$coefficients) <- c(min.tme,knots,max.tme)
  }
  rownames(return.value$coefficients) <- c("mu",paste("v",ind.unique,sep=""))
  return.value$fitted <- unlist(mapply(Xi,split(res.em$v,rep(1:res.em$n,each=res.em$p)),FUN=function(Xi,vi) Xi %*% (res.em$mu + vi),SIMPLIFY=FALSE))
  return.value$residuals <- y - as.vector(return.value$fitted)
  return.value$data <- data.frame(y=y,tme=tme,ind=ind)
  #return.value$em <- res.em
  return.value$call <- call
  return.value$dfMu <- res.em$dfMu
  return.value$dfV <- res.em$dfV
  return.value$logLik <- res.em$likelihood
  return.value$nobs <- length(y)
  return.value$df <- c(mu=res.em$dfMu,v=res.em$dfV)
  return.value$smoothingParameters <- c(mu=res.em$lambdaMu,v=res.em$lambdaV)
  return.value$parameters <- list(sigmaSquared=res.em$sigmaSquared,D=matrix(res.em$D,nrow=ncol(X)))
  return.value$iterations <- res.em$iterations
  return.value$info <- res.em$info
  if(!is.null(knots)) return.value$knots <- knots

  class(return.value) <- "sme"

  return(return.value)
}

sme.data.frame <- function(object,tme,ind,verbose=F,lambda.mu=NULL,lambda.v=NULL,maxIter=500,knots=NULL,zeroIntercept=F,deltaEM=1e-3,deltaNM=1e-3,criteria="AICc",...)
{
  if("variable" %in% names(object))
  {
    ys <- split(object$y,object$variable)
    tmes <- split(object$tme, object$variable)
    inds <- split(object$ind, object$variable)
    return(sme.list(ys,tmes,inds,verbose=verbose,lambda.mu=lambda.mu,lambda.v=lambda.v,maxIter=maxIter,knots=knots,zeroIntercept=zeroIntercept,deltaEM=deltaEM,deltaNM=deltaNM,criteria=criteria,...))
  }
  else
  {
    return(sme(object=object$y,tme=object$tme,ind=object$ind,,verbose=verbose,lambda.mu=lambda.mu,lambda.v=lambda.v,maxIter=maxIter,knots=knots,zeroIntercept=zeroIntercept,deltaEM=deltaEM,deltaNM=deltaNM,criteria=criteria,...))
  }
}

sme.list <- function(
  object,
  tme,
  ind,
  verbose=F,
  lambda.mu=NULL,
  lambda.v=NULL,
  maxIter=500,
  knots=NULL,
  zeroIntercept=F,
  deltaEM=1e-3,
  deltaNM=1e-3,
  criteria="AICc",
  numberOfThreads=-1,
  ...)
{
  ys <- object
  tmes <- tme
  inds <- ind

  for(i in 1:length(ys))
  {
    inds[[i]] <- as.factor(inds[[i]])
    ys[[i]] <- ys[[i]][order(inds[[i]])]
    tmes[[i]] <- tmes[[i]][order(inds[[i]])]
    inds[[i]] <- inds[[i]][order(inds[[i]])]
  }

  call <- match.call()

  if(!is.null(knots))
  {
    require(splines)
  
    max.tmes <- sapply(tmes,max)
    min.tmes <- sapply(tmes,min)
    allKnots <- mapply(max.tmes,min.tmes,FUN=function(max.tme,min.tme) c(min.tme,knots,max.tme),SIMPLIFY=FALSE)
    Bs <- lapply(allKnots,ns,knots=knots,intercept=T)

    Gs <- lapply(allKnots,function(a) roughnessMatrix(incidenceMatrix(a)))
    Gs <- mapply(Bs,Gs,FUN=function(B,G) t(B) %*% G %*% B,SIMPLIFY=FALSE)
  }
  else
  {
    Bs <- rep(NA,length(ys))
    allKnots <- rep(NA,length(ys))
    Gs <- lapply(tmes,function(tme) roughnessMatrix(incidenceMatrix(tme)))
  }

  M <- length(ys)
  Ni <- sapply(ys,length)

  yis <- mapply(ys,inds,FUN=function(y,ind) split(y,ind),SIMPLIFY=FALSE)
  ni <- sapply(yis, length)

  tmeis <- mapply(tmes,inds,FUN=function(tme,ind) split(tme,ind),SIMPLIFY=FALSE)
  Nij <- lapply(yis,sapply,length)

  if(is.null(knots))
  {
    Xs <- lapply(tmes, function(tme) incidenceMatrix(tme))
    Xis <- mapply(Xs,inds,FUN=function(X,ind) split.data.frame(X,ind),SIMPLIFY=FALSE)
  }
  else
  {
    Xs <- lapply(tmes,function(tme) ns(tme,knots=knots,intercept=T))
    Xis <- mapply(Xs,inds,FUN=function(X,ind) split.data.frame(X,ind),SIMPLIFY=FALSE)
  }

  ps <- sapply(Xs,ncol)

  if(any(is.null(lambda.mu) | is.null(lambda.v)))
  {
    res.em <- .C("SMEOptimizationMultiple",
                 M=as.integer(M),
                 y=as.double(unlist(ys)),
                 tme=as.double(unlist(tmes)),
                 ind=as.integer(unlist(inds)),
                 X=as.double(unlist(Xs)),
                 Ni=as.integer(Ni),
                 ni=as.integer(ni),
                 Nij=as.integer(unlist(Nij)),
                 p=as.integer(ps),
                 lambdaMu=double(M),
                 lambdaV=double(M),
                 G=as.double(unlist(Gs)),
                 mu=double(sum(ps)),
                 sigmaSquared=double(M),
                 D=double(sum(sapply(ps,"^",2))),
                 v=double(sum(mapply(ps,ni,FUN="*"))),
                 likelihood=double(M),
                 dfMu=double(M),
                 dfV=double(M),
                 zeroIntercept=as.integer(if(length(zeroIntercept)==1){ rep(zeroIntercept,M) }else{ zeroIntercept}),
                 iterations=integer(M),
                 maxIterations=as.integer(maxIter),
                 deltaEM=as.double(deltaEM),
                 deltaNM=as.double(deltaNM),
                 criteria=as.integer(match(criteria,c("AIC","AICc","BICN","BICn"))),
                 verbose=as.integer(verbose),
                 info=integer(M),
                 numberOfThreads=as.integer(numberOfThreads),
                 PACKAGE="sme")
  }
  else
  {
    if(length(lambda.mu)==1)
    {
      lambda.mu <- rep(lambda.mu,M)
    }
    if(length(lambda.v)==1)
    {
      lambda.v <- rep(lambda.v,M)
    }
    res.em <- .C("SMEMultiple",
                 M=as.integer(M),
                 y=as.double(unlist(ys)),
                 tme=as.double(unlist(tmes)),
                 ind=as.integer(unlist(inds)),
                 X=as.double(unlist(Xs)),
                 Ni=as.integer(Ni),
                 ni=as.integer(ni),
                 Nij=as.integer(unlist(Nij)),
                 p=as.integer(ps),
                 lambdaMu=as.double(lambda.mu),
                 lambdaV=as.double(lambda.v),
                 G=as.double(unlist(Gs)),
                 mu=double(sum(ps)),
                 sigmaSquared=double(M),
                 D=double(sum(sapply(ps,"^",2))),
                 v=double(sum(mapply(ps,ni,FUN="*"))),
                 likelihood=double(M),
                 dfMu=double(M),
                 dfV=double(M),
                 zeroIntercept=as.integer(if(length(zeroIntercept)==1){ rep(zeroIntercept,M) }else{ zeroIntercept}),
                 iterations=integer(M),
                 maxIterations=as.integer(maxIter),
                 deltaEM=as.double(deltaEM),
                 deltaNM=as.double(deltaNM),
                 verbose=as.integer(verbose),
                 info=integer(M),
                 numberOfThreads=as.integer(numberOfThreads),
                 PACKAGE="sme")

  }
  
  mus <- split(res.em$mu,rep(1:M,ps))
  vs <- split(res.em$v,rep(1:M,mapply(ps,ni,FUN="*")))
  iterations <- res.em$iterations
  dfMus <- res.em$dfMu
  dfVs <- res.em$dfV
  likelihoods <- res.em$likelihood
  Ds <- split(res.em$D,rep(1:M,ps*ps))
  sigmas <- res.em$sigmaSquared
  lambdaMus <- res.em$lambdaMu
  lambdaVs <- res.em$lambdaV
  infos <- res.em$info

  return.value <- mapply(
    mus,
    vs,
    Xs,
    Xis,
    Bs,
    ni,
    ps,
    allKnots,
    ys,
    tmes,
    inds,
    iterations,
    dfMus,
    dfVs,
    likelihoods,
    Ds,
    sigmas,
    lambdaMus,
    lambdaVs,
    infos,
    FUN=function(mu,v,X,Xi,B,n,p,allKnots,y,tme,ind,iterations,dfMu,dfV,likelihood,D,sigma,lambdaMu,lambdaV,info)
    {
      if(is.null(knots))
      {
        coefficients <- rbind(mu,t(matrix(v,ncol=n)))
        colnames(coefficients) <- attr(X,"tau")
      }
      else
      {
        coefficients <- rbind(as.vector(B %*% mu),do.call(rbind,lapply(split(v,rep(1:n,each=p)),function(v) as.vector(B %*% v))))
        colnames(coefficients) <- allKnots
      }
      rownames(coefficients) <- c("mu",paste("v",sort(unique(ind)),sep=""))
      fitted <- unlist(mapply(Xi,split(v,rep(1:n,each=p)),FUN=function(Xi,vi) Xi %*% (mu + vi),SIMPLIFY=FALSE))
      residuals <- y - as.vector(fitted)
      data <- data.frame(y=y,tme=tme,ind=ind)
      em <- list(iterations=iterations,likelihood=likelihood,dfMu=dfMu,dfV=dfV,D=matrix(D,nrow=p),sigmaSquared=sigma,lambdaMu=lambdaMu,lambdaV=lambdaV,info=info)
      return.value <- list(
        coefficients=coefficients,
        fitted=fitted,
        resid=residuals,
        data=data,
        #em=em,
        logLik=likelihood,
        df=c(mu=dfMu,v=dfV),
        call=call,
        nobs=length(y),
        smoothingParameters=c(mu=lambdaMu,v=lambdaV),
        parameters=list(sigmaSquared=sigma,D=matrix(D,nrow=p)),
        iterations=iterations,
        info=info)
      if(!is.null(knots)) return.value$knots <- knots
      class(return.value) <- "sme"
      return(return.value)
    },SIMPLIFY=FALSE)

  class(return.value) <- c("sme.list","list")
  return(return.value)
}

getRoughnessMatrix <- function(object)
{
  require(splines)

  if(is.null(object$knots))
  {
    X <- incidenceMatrix(object$data$tme)
    Xi <- split.data.frame(X,object$data$ind)
    G <- roughnessMatrix(X)
  }
  else
  {
    X <- ns(object$data$tme,knots=object$knots,intercept=T)
    Xi <- split.data.frame(X,object$data$ind)
  
    max.tme <- max(object$data$tme)
    min.tme <- min(object$data$tme)
    B <- ns(c(min.tme,object$knots,max.tme),knots=object$knots,intercept=T)
    G <- roughnessMatrix(incidenceMatrix(c(min.tme,object$knots,max.tme)))
    G <- t(B) %*% G %*% B
  }

  G
}

vcov.sme <- function(object,...)
{
  require(splines)

  if(is.null(object$knots))
  {
    X <- incidenceMatrix(object$data$tme)
    Xi <- split.data.frame(X,object$data$ind)
  }
  else
  {
    X <- ns(object$data$tme,knots=object$knots,intercept=T)
    Xi <- split.data.frame(X,object$data$ind)
  }

  G <- getRoughnessMatrix(object)
  Dv <- solve(solve(object$parameters$D) + object$smoothingParameters["v"] * G)
  Vi <- lapply(Xi,function(Xi) Xi %*% Dv %*% t(Xi) + diag(object$parameters$sigmaSquared,nrow=nrow(Xi)))
  inverseVi <- lapply(Vi,solve)
  inverseV <- blockDiagMat(inverseVi)
  tXinverseVX <- t(X) %*% inverseV %*% X
  conditionalCovarianceY <- solve(tXinverseVX + object$smoothingParameters["mu"] * G)
  
  vcov <- conditionalCovarianceY %*% tXinverseVX %*% conditionalCovarianceY
  rownames(vcov) <- colnames(coef(object))
  colnames(vcov) <- rownames(vcov)
  return(vcov)
}

rstandard.sme <- function(model,...)
{
  resid(model) / sqrt(model$parameters$sigmaSquared)
}

logLik.sme <- function(object,...)
{
  logLikelihood <- object$logLik
  attr(logLikelihood,"df") <- unname(object$df[1] + object$df[2])
  attr(logLikelihood,"nobs") <- length(resid(object))
  logLikelihood
}

AICc <- function(object)
{
  logLikelihood <- logLik(object)
  df <- attr(logLikelihood,"df")
  nobs <- nobs(object)
  AIC(object) + 2 * df * (df + 1) / (nobs - df - 1)
}

BICn <- function(object,...)
{
  UseMethod("BICn")
}

BICn.sme <- function(object,...)
{
  n <- nrow(coef(object))-1
  logLikelihood <- logLik(object,...)
  -2 * logLikelihood + log(n) * attr(logLikelihood,"df")
}

plot.sme <- function(x,type="model",...)
{
  if(type=="model")
  {
    plotSmeModel(x,...)
  }
  else if(type=="raw")
  {
    plotSmeRaw(x,...)
  }
  else if(type=="diagnostic")
  {
    plotSmeDiagnostic(x,...)
  }
  else
  {
    stop(paste("invalid plot type '",type,"'",sep=""))
  }
}

plotSmeModel <- function(x,xlab="Time",ylab="Y",showIndividuals=T,showConfidenceBands=F,col.meanCurve="red",...)
{
  mu <- spline(x=as.numeric(colnames(x$coefficients)),y=x$coefficients[1,],n=500,method="natural")
  fs <- lapply(2:nrow(x$coefficients),function(i){ spline(x=as.numeric(colnames(x$coefficients)),y=x$coefficients[1,] + x$coefficients[i,],method="natural",n=500) })
  ylim <- range(x$data$y, mu$y, sapply(fs,"[[","y"))

  if(showConfidenceBands)
  {
    mu.variance <- diag(vcov(x))
    upper.band <- spline(x=as.numeric(colnames(x$coefficients)),y=x$coefficients[1,] + 1.96 * sqrt(mu.variance),method="natural",n=500)
    lower.band <- spline(x=as.numeric(colnames(x$coefficients)),y=x$coefficients[1,] - 1.96 * sqrt(mu.variance),method="natural",n=500)
    ylim <- range(ylim, upper.band$y, lower.band$y)
  }

  plot(x=x$data$tme,y=x$data$y,ylim=ylim,xlab=xlab,ylab=ylab,...)
  if(showIndividuals)
  {
    for(i in 1:length(fs)){ lines(fs[[i]],lty="dashed") }
  }
  lines(mu,lwd=2,col=col.meanCurve)

  if(showConfidenceBands)
  {
    col.meanCurve.rgb <- col2rgb(col.meanCurve)
    polygon(x=c(upper.band$x,rev(lower.band$x)),y=c(upper.band$y,rev(lower.band$y)),col=rgb(col.meanCurve.rgb[1],col.meanCurve.rgb[2],col.meanCurve.rgb[3],alpha=125,maxColorValue=255),border=NA)
  }
}

plotSmeRaw <- function(x,xlab="Time",ylab="Y",mainTitle="",showModelFits=TRUE,showRawLines=FALSE)
{
  require(lattice)

  n <- nrow(coef(x))-1

  xyplot(y ~ tme | ind,
         data=x$data,
         xlab=xlab,
         ylab=ylab,
         main=mainTitle,
         subscripts=T,
         panel=
          function(...)
          {
            panel.xyplot(...)
            if(showRawLines)
            {
              panel.lines(...)
            }
            if(showModelFits)
            {
              ind <- x$data$ind[list(...)$subscripts[1]]
              f <- spline(x=as.numeric(colnames(coef(x))),y=coef(x)[1,] + coef(x)[paste("v",ind,sep=""),],method="natural",n=100)
              panel.lines(x=f$x,y=f$y,lty="dashed")
            }
          })
}

plotSmeDiagnostic <- function(x)
{
  layout(matrix(c(1,3,2,4),nrow=2))

  sres <- rstandard(x)

  plot(x=fitted(x),
       y=sres,
       xlab="Fitted values",
       ylab="Standardized residuals",
       main="Standardized residuals\nagainst fitted values")

  plot(x=x$data$tme,
       y=sres,
       xlab="Time",
       ylab="Standardized residuals",
       main="Standardized residuals\nagainst time")

  plot(x=x$data$y,
       y=sres,
       xlab="Response",
       ylab="Standardized residuals",
       main="Standardized residuals\nagainst response")

  qqnorm(sres,main="Q-Q plot for\nstandardized residuals")
  qqline(sres)
  layout(1)
}

blockDiagMat <- function(x)
{
     if(!is.list(x)) stop("x not a list")

     x <- x[as.logical(sapply(x, length))]

     n <- length(x)
     if(n==0) return(matrix(nrow=0,ncol=0))

     d <- array(unlist(lapply(x, dim)), c(2, n))
     rr <- d[1,]
     cc <- d[2,]
     rsum <- sum(rr)
     csum <- sum(cc)
     out <- array(0, c(rsum, csum))
     ind <- array(0, c(4, n))
     rcum <- cumsum(rr)
     ccum <- cumsum(cc)
     ind[1,-1] <- rcum[-n]
     ind[2,] <- rcum
     ind[3,-1] <- ccum[-n]
     ind[4,] <- ccum
     imat <- array(1:(rsum * csum), c(rsum, csum))
     iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
(y[3]+1):y[4]], imat=imat)
     iuse <- as.vector(unlist(iuse))
     out[iuse] <- unlist(x)
     return(out)
}

incidenceMatrix <- function(x)
{
  tau <- unique(x)
  tau <- sort(tau)
  M <- length(tau)

  X <- t(sapply(x,function(t){ as.numeric(t == tau) }))
 
  attr(X,"tau") <- tau
  attr(X,"M") <- M

  X
}

roughnessMatrix <- function(x)
{
  tau <- attr(x,"tau")
  K <- attr(x,"M")

  h <- vector(mode="numeric",length=K-1)

  for(r in 1:(K-1))
  {
    h[r] <- tau[r+1] - tau[r] 
  }

  A <- matrix(0, nrow=K, ncol=K-2)

  for(r in 1:(K-2))
  {
    A[r,r] <- 1 / h[r]

    A[r+1,r] <- -( (1 / h[r]) + (1 / h[r+1]) )

    A[r+2,r] <- 1 / h[r+1]
  }

  B <- matrix(0, nrow=K-2, ncol=K-2)

  B[1,1] <- (h[1] + h[2]) / 3

  B[2,1] <- h[2] / 6

  if(K > 4)
  {
    for(r in 1:(K-4))
    {
      B[r,r+1] <- h[r+1] / 6

      B[r+1,r+1] <- (h[r+1] + h[r+2]) / 3

      B[r+2,r+1] <- h[r+2] / 6
    }
  }

  B[K-3,K-2] <- h[K-2] / 6

  B[K-2,K-2] <- (h[K-2] + h[K-1]) / 3

  G <- A %*% solve(B) %*% t(A)
  attr(G,"A") <- A
  attr(G,"B") <- B
  G
}