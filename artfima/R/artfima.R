artfima <-
function(z, glp=c("ARTFIMA", "ARFIMA", "ARIMA"), arimaOrder=c(0,0,0), 
         constant=TRUE, likAlg=c("Whittle","exact"), 
         blueQ=FALSE, fixd=NULL, 
         optimMethod=c("LBFGSB","CG","NelderMead")) {
  #option fixd!=NULL only for ARTFIMA
  #
  stopifnot(is.numeric(z) && (is.ts(z) || is.vector(z)))
  stopifnot(length(arimaOrder==3) && is.numeric(arimaOrder) 
            && all(arimaOrder>=0))
  glp <- match.arg(glp)
  likAlg <- match.arg(likAlg)
  optimMethod <- match.arg(optimMethod)
  p <- arimaOrder[1]
  d0 <- arimaOrder[2]
  q <- arimaOrder[3]
  glpOrder <-switch(glp, "ARTFIMA"=2, "ARFIMA"=1, "ARIMA"=0)
  #fixd: must be null or numeric <2 and >=-0.5
  stopifnot(is.null(fixd)||(is.numeric(fixd)&&fixd<=2&&fixd>=-0.5)) 
  stopifnot(all(arimaOrder>=0))
  stopifnot(!(is.numeric(fixd)&&glpOrder!=2))
  #number of additional parameters, 0 for ARMA, 1 for ARFIMA, 2 for ARTFIMA
  #   except when fixd is not NULL then it is 1
  glpAdd <- glpOrder-ifelse(is.null(fixd), 0, 1)
  is.wholenumber <- function(x) abs(x - round(x)) < .Machine$double.eps^0.5
  stopifnot(is.wholenumber(p), is.wholenumber(d0), is.wholenumber(q))
#
  w <- z
  mnw <- mean(w)
  if(d0 > 0) {
    w <- diff(z, differences=d0)
    }
  if(!constant) {
    mnw <- 0
  }
  w <- w-mnw
  varw <- var(w)
  n <- length(w)
#initialization
  nbeta <- p+q+glpAdd #adjusted for fixd
  binit <- numeric(nbeta)
#Whittle method is fast because this is only done once
  if (likAlg=="Whittle") {
    Ip <- (spec.pgram(w, fast=FALSE, detrend=FALSE, plot=FALSE, taper=0)$spec)/(2*pi)
  }
  nullModelLoglikelihood <- (-n/2)*log(sum(w^2)/n)
  entropyPenalty <-  switch(likAlg,
      exact=-nullModelLoglikelihood,
      Whittle=sum(w^2)
  )
  entropyPenalty <- entropyPenalty+2*abs(entropyPenalty)
#negative log-likelihood=entropy
#begin Entropy. beta - pacf parameterization
#note - ***parameters passed using lexical scoping***
#also tacvf r is exported to arfima environment
  count <- 0
  Entropy<-function(beta) {
#in the optimization, put lambda, d, phi, theta
    phi<-theta<-lambda <- d <- numeric(0)
    count <<- count+1
    if (glpOrder==2) {
      lambda <- beta[1]
      d <- ifelse(is.null(fixd), beta[2], fixd) #fixd
    } else {
      if (glpOrder==1) d <- beta[1]
    }
    if(p>0) phi <- PacfToAR(beta[(1+glpAdd):(p+glpAdd)])
    if(q>0) theta <- PacfToAR(beta[(p+glpAdd+1):(p+q+glpAdd)]) 
    if (likAlg=="exact") {
      r <- tacvfARTFIMA(d=d, lambda=lambda, phi = phi, theta = theta, 
                        maxlag = n-1) 
#needed in some cases, eg. NileMin with p=1, q=3, glp="FGN"
      if (any(is.na(r))) {
        negLL <- NA
      } else {
        negLL <- try(-DLLoglikelihood(r, w), silent=TRUE)
      }
      negLL <- ifelse(is.numeric(negLL), negLL, entropyPenalty)
    } else {
      fp <- sdfartfima(n=n, d=d, lambda=lambda, phi=phi, theta=theta)
      negLL <- 2*mean(Ip/fp)
    }
#
#debugging...
#cat("\n ***********iter = ", iter, fill=TRUE)
#cat("count = ", count, fill=TRUE)
#cat("negLL=", negLL, fill=TRUE)
#cat("beta=",beta,fill=TRUE)
#cat("r[1:10] = ", r[1:10], fill=TRUE)
#cat("w[1:10] = ", w[1:10], fill=TRUE)
    negLL
  }#end Entropy()
#
#lower and upper limits with "L-BFGS-B"
  lambdaLo <- 0.0001
  lambdaHi <- 10
  dHi <- 10 #ARTFIMA limit
  dfHi <- 0.49 #ARFIMA limit
  if (glp=="ARTFIMA") {
    blo<- c(lambdaLo, -dHi, rep(-0.99,p+q))
    bhi<- c(lambdaHi, dHi,  rep(0.99,p+q))
  } else {
      if (glp=="ARFIMA") {
        blo<- c(-dfHi, rep(-0.99,p+q))
      } else { #ARMA
          blo<- rep(-0.99,p+q) 
        }
      bhi <- -blo
  }
#trace=6 for full output
  trace <- 0
#Brent when only 1 parameter
  if (length(binit)==1 && is.numeric(fixd)) {
    ans<-optim(par=binit, fn=Entropy, method="Brent", upper=lambdaHi, 
               lower=lambdaLo, control=list(trace=trace), hessian=TRUE)
    Alg <- "Brent"
  } 
    else {#GiantELSE
      Alg <- optimMethod
      ans <- switch(optimMethod, 
        CG = optim(par=binit, fn=Entropy, method="CG",
                    control=list(trace=trace, maxit=500), hessian=TRUE),
        LBFGSB = optim(par=binit, fn=Entropy, method="L-BFGS-B",
                 lower=blo, upper=bhi, control=list(trace=trace, maxit=500), 
                 hessian=TRUE),
        NelderMead = optim(par=binit, fn=Entropy, method="Nelder-Mead",
                       control=list(trace=trace, maxit=500), 
                       hessian=TRUE)
        )#end switch 
      }#endGiantELSE
      negLL <- ans$value
      bHat <- ans$par
      lambdaHat <- dHat <- phiHat <- thetaHat <- numeric(0)
      onBoundary <- FALSE
      if (glpOrder > 0) dHat <- bHat[glpOrder]
      if (glpOrder == 1 && abs(dHat) > dfHi) onBoundary <- TRUE
      if (glpOrder==2) {
        dHat <- ifelse(is.null(fixd), bHat[glpOrder], fixd)
        lambdaHat <-bHat[1]
        if (lambdaHat>=lambdaHi || lambdaHat <= lambdaLo) onBoundary <- TRUE
        if (abs(dHat) >= dHi) onBoundary <- TRUE
      }
    if (p > 0) phiHat <- PacfToAR(bHat[(1+glpAdd):(p+glpAdd)])
    if (q > 0) thetaHat <- PacfToAR(bHat[(p+1+glpAdd):(p+q+glpAdd)])
    rHat <- tacvfARTFIMA(d=dHat, lambda=lambdaHat, 
                         phi = phiHat, theta = thetaHat, maxlag = n-1)
    meanMLE <- 0
    if (blueQ) {
      meanMLE <- TrenchMean(rHat, w)
      w <- w-meanMLE
    }
  convergence <- ans$convergence
#end while###########################################################################  
#since Entropy is used, Hessian is pd
  Hinv<-try(solve(ans$hessian), silent=TRUE)
  if(!all(is.numeric(Hinv))) Hinv <- matrix(NA, nrow=nbeta, ncol=nbeta)
  if(!any(is.na(Hinv))){ #next square-root diagnonal elements
    sebHat <- suppressWarnings(sqrt(diag(Hinv)))
    } else {
      sebHat<-rep(NA, nbeta)
    }
  if (is.numeric(fixd)) sebHat <- c(sebHat[1],0,sebHat[-1])    
  rhoHat <- rHat[-1]/rHat[1]
  seMean <- sqrt(var(w)/n*(1+2*sum(1-(1:(n-1))/n*rhoHat^2))/n)
  constantHat <- mnw+meanMLE
  ansEx <- try(exactLoglikelihood(rHat, w), silent=TRUE)
  if(class(ansEx)!="try-error") { #possible converged but still problem
    LL <- ansEx$LL
    sigmaSq <- ansEx$sigmaSq
  } else {
    LL <- NA
    sigmaSq <- innovationVariance(w)
  }
  if (is.numeric(sebHat)&&likAlg=="Whittle") {
    sebHat <- sqrt(sigmaSq)*sebHat/sqrt(n)
  }
#
  snr <- (varw-sigmaSq)/sigmaSq
  res <- DLResiduals(rHat, w)
  out<-list(dHat=dHat, lambdaHat=lambdaHat, phiHat=phiHat, thetaHat=thetaHat, 
            constant=constantHat, seMean=seMean, se=sebHat, n=n, 
            sigmaSq=sigmaSq, snr=snr, likAlg=likAlg, convergence=convergence, 
            blueQ=blueQ, LL=LL, constantQ=constant, glp=glp, 
            arimaOrder=arimaOrder, glpOrder=glpOrder, fixd=fixd, glpAdd=glpAdd,
            tacvf=rHat, res=res, nullModelLogLik=nullModelLoglikelihood, 
            algorithm=Alg, 
            onBoundary=onBoundary, message=ans$message)
 class(out) <- "artfima"
 out
}
