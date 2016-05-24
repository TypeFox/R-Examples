profilemix.metaplus <- function(yi,sei,mods=NULL,justfit=FALSE,plotci=FALSE,slab=NULL) {
  
  isreg <- !is.null(mods)
  
  if (isreg) mods <- as.matrix(mods)
  
  ll.profilemix <- function(par,yi,sei,mods) {
    isreg <- !missing(mods)
    muhat <- par[1]
    tau2 <- par[2]
    tau2out <- par[3]
    lpoutlier <- par[4]
    if (isreg) xcoef <- matrix(par[5:length(par)],ncol=1)
    poutlier <- exp(lpoutlier)/(1+exp(lpoutlier))
    w <- 1.0/(tau2+sei^2)
    if (isreg) ll1 <- -0.5*(log(2*pi)-log(w)+w*(yi-muhat-as.vector(mods %*% xcoef))^2)+log(1-poutlier)
    else ll1 <- -0.5*(log(2*pi)-log(w)+w*(yi-muhat)^2)+log(1-poutlier)
    w <- 1.0/(tau2out+sei^2)
    if (isreg) ll2 <- -0.5*(log(2*pi)-log(w)+w*(yi-muhat-as.vector(mods %*% xcoef))^2)+log(poutlier)
    else ll2 <- -0.5*(log(2*pi)-log(w)+w*(yi-muhat)^2)+log(poutlier)
    l <- exp(cbind(ll1,ll2))
    negll <- -sum(log(apply(l,1,sum)))
    if (is.nan(negll)) negll <- NA
    if (!is.finite(negll)) negll <- NA
    if (is.na(negll)) negll <- 1e100
    return(negll)
  }      
  
  # obtain starting values
  if (isreg) {
    start.meta <- makestart.profilemix.metaplus(yi,sei, mods)$params
    start.val <- c(start.meta$muhat,start.meta$tau2,start.meta$tau2out,start.meta$lpoutlier,start.meta$xcoef)
    lower.val <- c(-Inf,0.0,0.0,-Inf,rep(-Inf,dim(mods)[2]))
  } else {
    start.meta <- makestart.profilemix.metaplus(yi,sei)$params
    start.val <- c(start.meta$muhat,start.meta$tau2,start.meta$tau2out,start.meta$lpoutlier)
    lower.val <- c(-Inf,0.0,0.0,-Inf)
  }     
  thenames <- c("muhat","tau2","tau2out","lpoutlier")
  if (isreg)  thenames <- c(thenames,dimnames(mods)[[2]])
  names(start.val) <- thenames
  names(lower.val) <- names(start.val)
  parnames(ll.profilemix) <- names(start.val)
  if (isreg) profilemix.fit <- mymle(ll.profilemix,start=start.val,vecpar=TRUE,optimizer="user",optimfun=myoptim,
                                     skip.hessian=TRUE,
                                     data=list(yi=yi,sei=sei,mods=mods),
                                     lower=lower.val)
  else profilemix.fit <- mymle(ll.profilemix,start=start.val,vecpar=TRUE,optimizer="user",optimfun=myoptim,
                               skip.hessian=TRUE,
                               data=list(yi=yi,sei=sei),
                               lower=lower.val)
  #print(profilemix.fit)
  results <- profilemix.fit@coef
  
  profilemix.profiled <- NULL
  
  if (!justfit)  {
    notprofiled <- TRUE
    while (notprofiled) {
      if (isreg) thehessian <- hessian(ll.profilemix,results,yi=yi,sei=sei,mods=mods)
      else thehessian <- hessian(ll.profilemix,results,yi=yi,sei=sei)
      # tau2 and tau2out can be zero 
      #browser()
      isproblem <- (results<1.0e-6) & ((1:length(results) %in% c(2,3)))
      isproblem2 <- isproblem*(1:length(results))
      noproblem2 <- (1-isproblem)*(1:length(results))
      if (all(isproblem2==0)) thehessian2 <- thehessian
      else thehessian2 <- thehessian[-isproblem2,-isproblem2]
      themyse <- suppressWarnings(sqrt(diag(ginv(thehessian2))))
      # expand back to original length
      myse <- rep(0,length(results))
      myse[noproblem2] <- themyse
      if (isreg) whichp <- c(1,5:(4+dim(mods)[2]))
      else whichp <- 1
      myse[is.na(myse)] <- 1.0e-6
      profilemix.stderr <- myse
      profilemix.profiled <- profilemix.profile(profilemix.fit,which=whichp,std.err=profilemix.stderr)
      if (class(profilemix.profiled) == "profile.mymle") notprofiled <- FALSE
      else {
        #browser()
        thenames <- c("muhat","tau2","tau2out","lpoutlier")
        start.val <- profilemix.profiled@fullcoef
        if (isreg) {
          lower.val <- c(-Inf,0.0,0.0,-Inf,rep(-Inf,dim(mods)[2]))
          thenames <- c(thenames,dimnames(mods)[[2]])
        } else {
          lower.val <- c(-Inf,0.0,0.0,-Inf)
        }
        #browser()
        parnames(ll.profilemix) <- thenames
        names(start.val) <- thenames
        names(lower.val) <- thenames
        if (isreg) profilemix.fit <- mymle(ll.profilemix,start=start.val,vecpar=TRUE,optimizer="user",
                                           data=list(yi=yi,sei=sei,mods=mods),
                                           skip.hessian=TRUE,
                                           control=list(eval.max=1000),
                                           lower=lower.val,optimfun=myoptim)
        else profilemix.fit <- mymle(ll.profilemix,start=start.val,vecpar=TRUE,optimizer="user",
                                     data=list(yi=yi,sei=sei),
                                     skip.hessian=TRUE,
                                     control=list(eval.max=1000),
                                     lower=lower.val,optimfun=myoptim)
        results <- profilemix.fit@coef
      }
    }

    if (any(order(profilemix.profiled@profile$muhat$z)!=(1:length(profilemix.profiled@profile$muhat$z)))) 
      warning("Profile loglikelihood is not unimodal in region of estimate. Possibly incorrect confidence intervals.")

    profilemix.ci <- confint(profilemix.profiled,method="uniroot")
    
    if (plotci) {
      tryCatch(plot(profilemix.profiled),
               error= function(e) {
                 #browser()
                 #plot(profilemix.profile@profile$muhat$z,profilemix.profile@profile$muhat$par.vals[,1])
                 print(paste("Error in CI plot: ",e))
               })
    }
    
    
    theci <- matrix(rep(NA,length(profilemix.fit@coef)*2),ncol=2)
    theci[whichp,] <- profilemix.ci
    results <- cbind(results,theci)
    
    pvalues <- rep(NA,length(start.val))
    for (iparm in whichp) {
      fixedparm <- names(start.val)[iparm]
      if (isreg) dostart <- paste("newstart.meta <- makestart.profilemix.metaplus(yi,sei, mods,",
                                  "fixed=list(",fixedparm,"=0.0))$params","\n",sep="")
      else dostart <- paste("newstart.meta <- makestart.profilemix.metaplus(yi,sei, NULL,",
                            "fixed=list(",fixedparm,"=0.0))$params","\n",sep="")
      eval(parse(text=dostart))
      if (isreg) newstart.val <- c(newstart.meta$muhat,newstart.meta$tau2,newstart.meta$tau2out,newstart.meta$lpoutlier,newstart.meta$xcoef)
      else newstart.val <- c(newstart.meta$muhat,newstart.meta$tau2,newstart.meta$tau2out,newstart.meta$lpoutlier)
      names(newstart.val) <- names(start.val)
      #      newstart.val <- newstart.val[-iparm]
      newlower.val <- lower.val
      if (isreg) doprofile <- paste("profilemix.fit0 <- mymle(ll.profilemix,start=newstart.val,vecpar=TRUE,\n",
                                    "optimizer=\"user\",optimfun=myoptim,data=list(yi=yi,sei=sei,mods=mods),\n",
                                    "skip.hessian=TRUE,\n",
                                    "lower=newlower.val,fixed=list(",fixedparm,"=0.0))",sep="")
      else doprofile <- paste("profilemix.fit0 <- mymle(ll.profilemix,start=newstart.val,vecpar=TRUE,\n",
                              "optimizer=\"user\",optimfun=myoptim,data=list(yi=yi,sei=sei),\n",
                              "skip.hessian=TRUE,\n",
                              "lower=newlower.val,fixed=list(",fixedparm,"=0.0))",sep="")
      eval(parse(text=doprofile))
      pvalues[iparm] <- anova(profilemix.fit,profilemix.fit0)[2,5]
    }
    results <- cbind(results,pvalues)
    dimnames(results)[[2]] <- c("Est.","95% ci.lb","95% ci.ub","pvalue")
    # transform lpoutlier to poutlier
    dimnames(results)[[1]][4] <- "Outlier prob."
    results[4,1] <- 1.0/(1+exp(-results[4,1]))
  }

  # calculate final posterior probabilities
  muhat <- profilemix.fit@coef[1]
  tau2 <- profilemix.fit@coef[2]
  tau2out <- profilemix.fit@coef[3]
  lpoutlier <- profilemix.fit@coef[4]
  if (isreg) xcoef <- matrix(profilemix.fit@coef[5:length(profilemix.fit@coef)],ncol=1)
  
  poutlier <- exp(lpoutlier)/(1+exp(lpoutlier))
  
  w <- 1.0/(tau2+sei^2)
  if (isreg) ll1 <- -0.5*(log(2*pi)-log(w)+w*(yi-muhat-as.vector(mods %*% xcoef))^2)
  else ll1 <- -0.5*(log(2*pi)-log(w)+w*(yi-muhat)^2)
  w <- 1.0/(tau2out+sei^2)
  if (isreg) ll2 <- -0.5*(log(2*pi)-log(w)+w*(yi-muhat-as.vector(mods %*% xcoef))^2)
  else ll2 <- -0.5*(log(2*pi)-log(w)+w*(yi-muhat)^2)
  p <- c(1-poutlier,poutlier)
  l <- exp(cbind(ll1,ll2))
  
  prop <- t(p*t(l))/apply(t(p*t(l)),1,sum)
  
  return(list(results=results,yi=yi,sei=sei,mods=mods,slab=slab,fittedmodel=profilemix.fit,
              outlier.prob=prop[,2],justfit=justfit,profile=profilemix.profiled,random="mixture"))
}