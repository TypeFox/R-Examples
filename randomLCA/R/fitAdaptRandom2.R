fitAdaptRandom2 <- function(outcomes,freq,nclass=2,initoutcomep,initclassp,initlambdacoef,initltaucoef,
                                        level2size,constload,calcSE=FALSE,gh,probit,byclass,qniterations,penalty,verbose=FALSE) {
  
  #	print(initoutcomep)
  #	print(initclassp)
  #	print(initlambdacoef)
  #	print(initltaucoef)
  
  # parameters
  #   outcomes matrix of outcomes 0 or 1
  #   freq vector of frequencies corresponding to each outcome combination
  #   nclass number of classes
  #   initoutcomep initial outcome probabilities
  #   initclassp initial class probabilities
  #   initlambdacoef initial lambda coefficient
  #   initltaucoef initial tau coefficient log scale
  #   level2size number of outcomes in each period
  #   constload vary loading for each outcome
  #   calcSE calculate standard errors ?
  #   gh matrix of gauss-hermite coefficients first column positions, second columns weights
  #   probit use probit transform rather than logitic to obtain outcome probabilities
  #   verbose print information about algorithm    
  
  outcomes <- as.matrix(outcomes)
  mode(outcomes) <- "double"
  
  nlevel1 <- level2size
  nlevel2 <- dim(outcomes)[2]/level2size
  nlevel3 <- length(freq)
  
  if (constload) nlambda <- 1
  else nlambda <- nlevel1
  
  outcomestart <- nclass
  outcomeend <- (nclass+nlevel1*nlevel2*nclass-1)
  
  if (byclass) {
    lambdastart <- (nlevel1*nlevel2*nclass+nclass)
    lambdaend <- (nlevel1*nlevel1*nclass+nclass+nclass*nlambda-1)
    taustart <- (nlevel1*nlevel2*nclass+nclass+nclass*nlambda)
    tauend <- (nlevel1*nlevel2*nclass+nclass+nclass*nlambda+nclass-1)
  } else {
    lambdastart <- (nlevel1*nlevel2*nclass+nclass)
    lambdaend <- (nlevel1*nlevel2*nclass+nclass+nlambda-1)
    taustart <- (nlevel1*nlevel2*nclass+nclass+nlambda)
    tauend <- (nlevel1*nlevel2*nclass+nclass+nlambda)
  }
  
  calclikelihood <- function(classx,outcomex,lambdacoef,ltaucoef,momentdata,gh,
                             updatemoments=FALSE,calcfitted=FALSE,zprop=NULL) {
    
#       print("classx")
#     print(classx)
#     print("outcomex")
#     print(outcomex)
#       print("lambdacoef")
#       print(lambdacoef)
#       print("ltaucoef")
#       print(ltaucoef)
    
    # turn classp into actual probabilities
    classp2 <- c(0,classx)       
    classp2 <- exp(classp2)/sum(exp(classp2))
    
    #	print(outcomex)
    #	print(classp2)
    #	print(lambdacoef)
    #	print(exp(ltaucoef))
    
    #		browser()
    
    newmoments <- NULL
    ill <- matrix(rep(NA,nclass*nlevel3),ncol=nclass)
    mylambdacoef <- lambdacoef
    if (constload) {
      if (byclass) mylambdacoef <- matrix(rep(lambdacoef,nlevel1),nrow=nclass)
      else mylambdacoef <- rep(lambdacoef,nlevel1)       
    }
    #   print("mylambdacoef")
    #   print(mylambdacoef)
    for (iclass in 1:nclass) {
      result <- .Call("bernoulliprobrandom2",outcomes,outcomex[iclass,],
                      if (byclass) mylambdacoef[iclass,]  else mylambdacoef,
                      if (byclass) ltaucoef[iclass] else ltaucoef,
                      gh,momentdata[,((iclass-1)*(2+nlevel2*3)+1):(iclass*(2+nlevel2*3))],
                      probit,updatemoments)
      ill[,iclass] <- exp(result[[1]])*classp2[iclass]
      if (updatemoments) newmoments <- cbind(newmoments,result[[2]])
    }
    #		browser()
    # if zprop not supplied then we have the usual maximum likelihood
    if (is.null(zprop)) {
      ill2 <- log(rowSums(ill))
      ll <- sum(ill2*freq)
    } else {
      # otherwise calculate the comple data maximum likelihood for the em algorithm
      #        browser()
      ill2 <- rowSums(log(ill)*zprop)
      ll <- sum(ill2*freq)			  
    }
    # penalise extreme outcome probabilities
    outcomep <- as.vector(1/(1+exp(abs(outcomex))))
#    pen <- dbeta(outcomep,1+penalty,1+penalty,log=TRUE)
    pen <-SciencesPo::ddirichlet(matrix(outcomep,nrow=1),rep(1+penalty/(nclass*2),length(outcomep)),log=TRUE)-SciencesPo::ddirichlet(matrix(outcomep,nrow=1),rep(1,length(outcomep)),log=TRUE)
    penll <- ll+sum(pen)
    if (is.nan(penll) || is.infinite(penll)) penll <- -1.0*.Machine$double.xmax
    if (calcfitted) {
      fitted <- exp(ill2)*sum(ifelse(apply(outcomes,1,function(x) 
        any(is.na(x))),0,freq))*
        ifelse(apply(outcomes,1,function(x) any(is.na(x))),NA,1)
      classprob <- ill/exp(ill2)
      return(list(logLik=ll,penlogLik=penll,moments=newmoments,fitted=fitted,classprob=classprob))
    } else return(list(logLik=ll,penlogLik=penll,moments=newmoments))
  }  # end of calclikelihood
  
  adaptivefit <- function(classx,outcomex,lambdacoef,ltaucoef,momentdata,gh) {
    
    fitparams <- function(classx,outcomex,lambdacoef,ltaucoef,
                          momentdata,gh,calcSE,noiterations=qniterations,zprop=NULL) {
      
      calcllfornlm <- function(params,momentdata,gh,zprop) {  
        oneiteration <- calclikelihood(if (nclass==1) NULL else params[1:(nclass-1)],
                                       matrix(params[outcomestart:outcomeend],nrow=nclass),
                                       if (byclass) matrix(params[lambdastart:lambdaend],nrow=nclass) else params[lambdastart:lambdaend],
                                       params[taustart:tauend],
                                       momentdata,gh,zprop=zprop)
        return(-oneiteration$penlogLik)
      }
      
      nlm1 <- nlm(calcllfornlm, c(classx,as.vector(outcomex),lambdacoef,ltaucoef),
                  iterlim = noiterations,
                  print.level=ifelse(verbose,2,0),
                  check.analyticals = FALSE,hessian=calcSE,momentdata=momentdata,gh=gh,zprop=zprop)
      return(list(penlogLik=-nlm1$minimum,
                  classx=if (nclass==1) NULL else nlm1$estimate[1:(nclass-1)],
                  outcomex=matrix(nlm1$estimate[outcomestart:outcomeend],nrow=nclass),
                  lambdacoef= if (byclass) matrix(nlm1$estimate[lambdastart:lambdaend],nrow=nclass) else nlm1$estimate[lambdastart:lambdaend],
                  ltaucoef=nlm1$estimate[taustart:tauend],
                  nlm=nlm1))
    }
    
    
    oneiteration <- calclikelihood(classx, outcomex, lambdacoef,ltaucoef,
                                   momentdata,gh,updatemoments=TRUE)
    currll <- oneiteration$penlogLik
    if (verbose) cat('Initial ll',currll,"\n")
    lastll <- 2*currll
    # shift the quadrature points for the first time
    while (abs((lastll-currll)/lastll)>1.0e-6) {
      lastll <- currll
      momentdata <- oneiteration$moments
      oneiteration <- calclikelihood(classx,outcomex,lambdacoef,ltaucoef,
                                     momentdata,gh,updatemoments=TRUE)
      currll <- oneiteration$penlogLik
      zprop <- oneiteration$classprobs
      if (verbose) cat("current ll",currll,"\n")       
    }
    
    adaptive <- TRUE
    prevll <- -Inf
    nadaptive <- 0
    while(adaptive) {
      # need to do an optimisation on the other parameters
      fitresults <- fitparams(classx,outcomex,lambdacoef,ltaucoef,
                              momentdata,gh,calcSE=FALSE,zprop=zprop)
      currll <- fitresults$penlogLik
      outcomex <- fitresults$outcomex
      classx <- fitresults$classx
      lambdacoef <- fitresults$lambdacoef
      ltaucoef <- fitresults$ltaucoef
      if (verbose) cat("current ll from optimisation",currll,"\n") 
      optll <- currll
      # shift the quadrature points again
      oneiteration <- calclikelihood(classx,outcomex,lambdacoef,ltaucoef,
                                     momentdata,gh,updatemoments=TRUE)
      currll <- oneiteration$penlogLik
      lastll <- 2*currll
      while(abs((lastll-currll)/lastll)>1.0e-7) {
        lastll <- currll
        momentdata <- oneiteration$moments
        oneiteration <- calclikelihood(classx,outcomex,lambdacoef,ltaucoef,
                                       momentdata,gh,updatemoments=TRUE)
        currll <- oneiteration$penlogLik
        if (verbose) cat("current ll",currll,"\n")       
      }
      adaptive <- (abs((currll-optll)/currll)>1.0e-7) ||
        (abs((currll-prevll)/currll)>1.0e-7)
      if ((prevll-currll)/abs(currll) > 1.0e-3) stop("divergence - increase quadrature points")
      nadaptive <- nadaptive+1
      if (nadaptive > 200) stop("too many adaptive iterations - increase quadrature points")
      prevll <- currll
      # if (is.na(adaptive)) browser()
      zprop <- oneiteration$classprobs
    }
    fitresults <- fitparams(classx,outcomex,lambdacoef,ltaucoef,
                            momentdata,gh,calcSE=calcSE,noiterations=500)
    return(list(nlm=fitresults$nlm,momentdata=momentdata))
  } # end adaptivefit
  
  # momentdata is level3 and level2
  # mu3,tau3,nlevel2*(mu2,taucoef,gamma2)
  # repeated for each class
  
  momentdata <- matrix(rep(c(rep(c(0,1),each=nlevel3),
                             rep(rep(c(0,1,0),each=nlevel3),times=nlevel2)),times=nclass),
                       nrow=nlevel3)
  
  if (nclass==1) classx <- NULL
  else  {
    classx <- rep(NA,nclass-1)
    initclassp <- ifelse(initclassp<1.0e-4,1.0e-4,initclassp)        
    initclassp <- ifelse(initclassp>(1.0-1.0e-4),1-1.0e-4,initclassp)          
    for (i in 2:nclass) classx[i-1] <- log(initclassp[i]/initclassp[1])
  }
  
  initoutcomep <- ifelse(initoutcomep<1.0e-4,1.0e-4,initoutcomep)
  initoutcomep <- ifelse(initoutcomep>(1.0-1.0e-4),1-1.0e-4,initoutcomep)
  if (probit) outcomex <- qnorm(initoutcomep)
  else outcomex <- log(initoutcomep/(1-initoutcomep))
  
  if (missing(initlambdacoef) || is.null(initlambdacoef)) {
    if (byclass) lambdacoef <- matrix(rep(0,level2size*nclass),nrow=nclass)
    else lambdacoef <- rep(0,level2size)
  } else lambdacoef <- initlambdacoef
  
  # choose among possible ltaucoef
  if (missing(initltaucoef) || is.null(initltaucoef)) {
    testltaucoef <- -3.0
    maxltau <- NA
    maxll <- -Inf
    repeat {
      if (verbose) cat('trying ltaucoef ',testltaucoef,"\n")
      if (byclass) theltaucoef <- rep(testltaucoef,nclass)
      else theltaucoef <- testltaucoef
      # browser()
      onelikelihood <- calclikelihood(classx,outcomex,lambdacoef,theltaucoef,momentdata,gh)
      currll <- onelikelihood$penlogLik
      if (verbose) cat("ll",currll,"\n")		
      # when the ll starts decreasing, give up
      #       print("currll")
      #       print(currll)
      if (currll < maxll) break()
      maxll <- currll
      maxltau <- testltaucoef
      testltaucoef <- testltaucoef+0.1
    }
    if (verbose) cat('using ltaucoef ',maxltau,"\n")
    if (byclass) ltaucoef <- rep(maxltau,nclass)
    else ltaucoef <- maxltau
  }
  else ltaucoef <- initltaucoef
  
  myfit <- adaptivefit(classx, outcomex, lambdacoef,ltaucoef,momentdata,gh)
  
  optim.fit <- myfit$nlm
  momentdata <- myfit$momentdata
  
  # extract the se
  #browser()
  if (!calcSE) separ <- rep(NA,length(optim.fit$estimate))
  else {
    s <- svd(optim.fit$hessian)
    separ <- sqrt(diag(s$v %*% diag(1/s$d) %*% t(s$u)))
    separ[!is.finite(separ)] <- NA
  }
  # calculate the probabilities
  if (nclass==1) classp <- 1
  else {
    classp <-optim.fit$estimate[1:(nclass-1)]
    # add extra column to classp
    classp <- c(0,classp)       
  }       
  outcomex <- matrix(optim.fit$estimate[outcomestart:outcomeend],nrow=nclass)
  # transform using logistic to probabilities     
  classp <- exp(classp)/sum(exp(classp))
  if (probit) outcomep <- pnorm(outcomex)
  else outcomep <- exp(outcomex)/(1+exp(outcomex))
  if (byclass) lambdacoef <- matrix(optim.fit$estimate[lambdastart:lambdaend],nrow=nclass)
  else lambdacoef <- optim.fit$estimate[lambdastart:lambdaend]
  ltaucoef <- optim.fit$estimate[taustart:tauend]
  
  calcrandom <- function() {
    
    outcomex <- log(outcomep/(1-outcomep))
    
    onerandom <- function(x) {
      
      loglik <- function(beta) {
        # calculate probabilities under each class
        for (i in 1:nclass) {
          # calculate the outcome probabilities for this class and current random
          if (byclass) {
            if (probit) outcomep <- pnorm(outcomex[i,]+rep(lambdacoef[i,],nlevel2)*
                                            (beta[1]+exp(ltaucoef[i])*rep(beta[2:(1+nlevel2)],each=nlevel1)))
            else  outcomep <- 1/(1+exp(-outcomex[i,]-rep(lambdacoef[i,],nlevel2)*
                                         (beta[1]+exp(ltaucoef[i])*rep(beta[2:(1+nlevel2)],each=nlevel1))))
          } else {
            if (probit) outcomep <- pnorm(outcomex[i,]+rep(lambdacoef,nlevel2)*
                                            (beta[1]+exp(ltaucoef)*rep(beta[2:(1+nlevel2)],each=nlevel1)))
            else outcomep <- 1/(1+exp(-outcomex[i,]-rep(lambdacoef,nlevel2)*
                                        (beta[1]+exp(ltaucoef)*rep(beta[2:(1+nlevel2)],each=nlevel1))))
          }
          oneprob <- t(apply(t(x)*outcomep+t(1-x)*(1-outcomep),2,prod,na.rm=TRUE))
          # multiply by class probabilities
          if (i==1) allprob <- oneprob*classp[i]
          else allprob <- allprob+oneprob*classp[i]
        }
        
        ll <- -(sum(log(allprob))+sum(dnorm(beta,mean=0,sd=1,log=TRUE)))
        if (is.nan(ll) || is.infinite(ll)) ll <- 1.0*.Machine$double.xmax
        return(ll)
      }
      beta <- rep(0,1+nlevel2)
      optim.fit <- nlm(loglik,beta,print.level=0,iterlim=1000,hessian=TRUE,gradtol=1.0e-7)
      beta <- optim.fit$estimate
      sebeta <- sqrt(diag(solve(optim.fit$hessian)))
      checkx <- matrix(x,ncol=nlevel1,byrow=T)
      checkx <- apply(checkx,1,function(x) any(is.na(x)))
      checkx <- c(FALSE,checkx)
      beta[checkx] <- NA
      sebeta[checkx] <- NA
      return(c(beta=beta,sebeta=sebeta))
    }
    betas <- t(apply(outcomes,1,onerandom))
    return(betas)
  }
  
  ranef <- calcrandom()
  
  if (byclass) final <- calclikelihood(if (nclass==1) NULL else optim.fit$estimate[1:(nclass-1)],
                                       matrix(optim.fit$estimate[outcomestart:outcomeend],nrow=nclass),
                                       matrix(optim.fit$estimate[lambdastart:lambdaend],nrow=nclass),
                                       optim.fit$estimate[taustart:tauend],
                                       momentdata,gh,updatemoments=FALSE,calcfitted=TRUE)    
  else final <- calclikelihood(if (nclass==1) NULL else optim.fit$estimate[1:(nclass-1)],
                               matrix(optim.fit$estimate[outcomestart:outcomeend],nrow=nclass),
                               optim.fit$estimate[lambdastart:lambdaend],
                               optim.fit$estimate[taustart:tauend],
                               momentdata,gh,updatemoments=FALSE,calcfitted=TRUE)    
  fitted <- final$fitted
  classprob <- final$classprob
  
  np <- length(optim.fit$estimate)
  nobs <- sum(freq)
  deviance <- 2*sum(ifelse(freq==0,0,freq*log(freq/fitted)))
  
  # check that results are sensible
  if (any(abs(as.vector(outcomex))>20)) warning("Problem in solution, probably unstable")
  
  list(fit=optim.fit,nclass=nclass,classp=classp,outcomep=outcomep,lambdacoef=lambdacoef,taucoef=exp(ltaucoef),se=separ,
       np=np,nobs=nobs,logLik=final$logLik,penlogLik=final$penlogLik,freq=freq,fitted=fitted,ranef=ranef
       ,classprob=classprob,deviance=deviance)
}
