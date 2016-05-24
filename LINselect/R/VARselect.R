# ---------------------------------------------------
VARselect <- function
  ### Estimation in the regression model : \eqn{Y= X \beta + \sigma N(0,1)}\cr
  ### Variable selection by choosing the best predictor among
  ### predictors emanating \cr from different methods as lasso,
  ### elastic-net, adaptive lasso, pls, randomForest.
  ##references<<  See Baraud  et al. 2010 
  ## \url{http://hal.archives-ouvertes.fr/hal-00502156/fr/} \cr
  ## Giraud et al., 2013,
  ## \url{http://projecteuclid.org/DPubS?service=UI&version=1.0&verb=Display&handle=euclid.ss/1356098553}
(Y, ##<<vector with n components : response variable. 
 X, ##<<matrix with n rows and p columns : covariates.
 dmax=NULL,##<<integer : maximum number of variables in the lasso
 ##estimator. \code{dmax} \eqn{\le} D where \cr 
 ##  D = min (3*p/4 , n-5) if  p\eqn{ \ge }n
 ## \cr D= min(p,n-5) if
 ## p < n. \cr Default : \code{dmax} = D.
 normalize=TRUE, ##<<logical : if TRUE the columns of X are scaled 
 method=c("lasso","ridge","pls","en","ALridge","ALpls","rF","exhaustive"), ##<< vector of characters  whose components are subset of \cr
  ## \dQuote{lasso}, \dQuote{ridge}, \dQuote{pls}, \dQuote{en},
 ## \dQuote{ALridge}, \dQuote{ALpls}, \dQuote{rF},
 ## \dQuote{exhaustive}.
 pen.crit=NULL, ##<<vector with \code{dmax}+1 components : for d=0,
 ##..., \code{dmax}, \code{penalty[d+1]} gives the value of the
 ##penalty for the dimension d. Default : \code{penalty} = NULL. In
 ##that case, the
 ##penalty will be calculated by the function penalty. 
 lasso.dmax=NULL, ##<<integer lower than \code{dmax}, default = \code{dmax}.
 ridge.dmax=NULL, ##<<integer lower than \code{dmax}, default = \code{dmax}.
 pls.dmax=NULL, ##<<integer lower than \code{dmax}, default = \code{dmax}.
 en.dmax=NULL, ##<<integer lower than \code{dmax}, default = \code{dmax}.
 ALridge.dmax=NULL, ##<<integer lower than \code{dmax}, default = \code{dmax}.
 ALpls.dmax=NULL, ##<<integer lower than \code{dmax}, default = \code{dmax}.
 rF.dmax=NULL, ##<<integer lower than \code{dmax}, default = \code{dmax}.
 exhaustive.maxdim=500000, ##<<integer : maximum number of subsets of
 ##covariates considered  in the exhaustive method. See details.
 exhaustive.dmax=NULL, ##<<integer lower than \code{dmax}, default = \code{dmax}
 en.lambda=c(0.01,0.1,0.5,1,2,5), ##<<vector : tuning parameter of the
 ##ridge. It is the input parameter \code{lambda} of function
 ##\code{\link[elasticnet]{enet}}
 ridge.lambda=c(0.01,0.1,0.5,1,2,5), ##<<vector : tuning parameter of the
 ##ridge. It is the input parameter lambda of function
 ##\code{\link{lm.ridge}}
 rF.lmtry=2, ##<<vector : tuning paramer \code{mtry} of function
 ##\code{\link[randomForest]{randomForest}}, \code{mtry} =p/\code{rF.lmtry}.
 pls.ncomp=5,  ##<<integer : tuning parameter of the pls. It is the
 ##   input parameter \code{ncomp} of the function
 ##   \code{\link[pls:mvr]{plsr}}. See details.
 ALridge.lambda=c(0.01,0.1,0.5,1,2,5), ##<< similar to
 ## \code{ridge.lambda} in the adaptive lasso procedure.
 ALpls.ncomp=5,  ##<< similar to \code{pls.ncomp} in the
 ##adaptive lasso procedure. See details.
 max.steps=NULL, ##<<integer. Maximum number of steps in the lasso
 ##procedure. Corresponds to the input \code{max.steps} of the function
 ##\code{\link[elasticnet]{enet}}. \cr
 ##Default :
 ##\code{max.steps} = 2*min(p,n)
 K=1.1,##<<scalar : value of the parameter \eqn{K} in the LINselect criteria.
 verbose=TRUE,##<<logical : if TRUE a trace of the current process is displayed in real time.
 long.output=FALSE##<<logical : if FALSE only the component summary
 ## will be returned. See Value.
 ) {
  #
  n <- length(Y)
  if (normalize) X <- scale(X)
  # penalty calculation
  result <- NULL
  p <- dim(X)[2]
  if (p>=n) dmax <- floor(min(c(3*p/4,n-5,dmax)))
  if (p<n) dmax <- floor(min(c(p,n-5,dmax)))
  if (!is.null(pen.crit)) {
    if (length(pen.crit)==(dmax+1)) pen=pen.crit
    if (length(pen.crit)!=(dmax+1)) {
      print(paste("warning: the length of pen.crit is not equal to",dmax+1))
      print(paste("the value of pen.crit is set to NULL"))
      pen.crit=NULL
    }
  }
   if (is.null(pen.crit)) {
    D <- 0:dmax
    dm <- pmin(D,rep(p/2,dmax+1))
    Delta <- lgamma(p+1)-lgamma(dm+1)-lgamma(p-dm+1)+2*log(D+1)
    pen <- penalty(Delta, n, p, K)
    result$pen.crit=pen
  }
  crit <- NULL
  supp <- list(NULL)
  meth <- NULL
  ll=0
  if (is.null(max.steps)) max.steps = min(n,p)
  # lasso procedure
  if (is.element("lasso",method)) {
    if (verbose==TRUE) print(
          paste("LASSO PROCEDURE with option max.steps=",max.steps,"dmax=",min(lasso.dmax, dmax)))
#    library(elasticnet)
    ##note<< When method is \code{lasso}, library \code{elasticnet} is loaded.
    res.lasso <- try(enet(X,Y,lambda=0,intercept=TRUE,normalize=FALSE,max.steps=max.steps))
    if (inherits(res.lasso, "try-error")) {
      print("error with method=lasso, when calling enet")
      result$lasso <- res.lasso[1]
    }
    if (!inherits(res.lasso, "try-error")) {
      BetaHat <- array(0,c(dim(res.lasso$beta.pure)[1],p))
      BetaHat[,res.lasso$allset] <- res.lasso$beta.pure
      DimBetaHat <- apply(BetaHat!=0,1,sum)
      Nmod.lasso  <- dim(BetaHat)[1]
      ff <- calc.fhat(Y,X,BetaHat,Nmod.lasso,n,p)
      f.proj <- ff$fhat
      r.proj <- ff$r
      m.lasso <- ff$m
      Nmod <-  which(DimBetaHat<=min(lasso.dmax, dmax))
      res1 <- EstSelect5(Y,X,f.proj[,Nmod],m.lasso[Nmod],r.proj[Nmod], pen)
      if (long.output) {
        result$chatty$lasso<-list(support=m.lasso[Nmod],crit=res1$crit)
      }
      result$lasso<-list(support=res1$mHat,crit=res1$crit[res1$lHat],fitted=res1$fHat)
      crit <- c(crit,res1$crit[res1$lHat])
      ll=ll+1
      supp[[ll]]<-res1$mHat
      meth <- c(meth,"lasso")
    }
  }
  # elastic net procedure
  if (is.element("en",method)) {
    if (verbose==TRUE) {
      print(paste(c("ELASTIC-NET PROEDURE with options max.steps=",max.steps,
            "dmax=",min(en.dmax, dmax),"en.lambda=",en.lambda),collapse=" "))
    }
#    library(elasticnet)
    ##note<< When method is \code{en}, library \code{elasticnet} is loaded.
    m.lasso <- list(NULL)
    Nm <- 0
    f.proj<-NULL
    r.proj <- NULL
    nberr <- 0
    for (iL in 1:length(en.lambda)) {
      res.lasso <-
        try(enet(X,Y,lambda=en.lambda[iL],intercept=TRUE,normalize=FALSE,max.steps=max.steps))
      if (inherits(res.lasso, "try-error")) {
        nberr=nberr+1
        print(paste("error with method=en, when calling enet with parameter
lambda =",en.lambda[iL]))
        next
      }
      BetaHat <- array(0,c(dim(res.lasso$beta.pure)[1],p))
      BetaHat[,res.lasso$allset] <- res.lasso$beta.pure
      DimBetaHat <- apply(BetaHat!=0,1,sum)
      Nmod.lasso  <- dim(BetaHat)[1]
      ff <- calc.fhat(Y,X,BetaHat,Nmod.lasso,n,p)
      f.projiL <- ff$fhat
      r.projiL <- ff$r
      m.lassoiL <- ff$m
      Nmod <-  which(DimBetaHat<=min(en.dmax, dmax))
      f.proj <- cbind(f.proj,f.projiL[,Nmod])
      r.proj <- c(r.proj,r.projiL[Nmod])
      m.lasso[(Nm+1):(Nm+length(Nmod))]  <- m.lassoiL[Nmod]
      Nm<- Nm+length(Nmod)
    }
    if (nberr==length(en.lambda)) {
      result$en <- "error in enet for all values of lambda"
    }
    if (nberr!=length(en.lambda)) {   
      res1 <-  EstSelect5(Y, X, f.proj, m.lasso, r.proj, pen)
      if (long.output) {
        result$chatty$en <-list(support=m.lasso,crit=res1$crit)
      }
      result$en<-list(support=res1$mHat,crit=res1$crit[res1$lHat],fitted=res1$fHat)
      crit <- c(crit,res1$crit[res1$lHat])
      ll=ll+1
      supp[[ll]]<-res1$mHat
      meth <- c(meth, "en")
    }
  }
  # ridge procedure
  if (is.element("ridge",method)) {
    if (verbose==TRUE) {
      print(paste(c("RIDGE PROCEDURE with options dmax=",
                    min(ridge.dmax,dmax),
                    "ridge.lambda=",ridge.lambda),collapse=" "))
    }
#     library(MASS)
    ##note<< When method is \code{ridge}, library \code{MASS} is loaded.
    dmax.rdg <- min(ridge.dmax,dmax)
    m.rdg <- list(NULL)
    f.proj<-NULL
    r.proj <- NULL
    res.rdg <- try(lm.ridge(Y~X,lambda=ridge.lambda))
    if (inherits(res.rdg, "try-error")) {
      print("error with method=ridge, when calling lm.ridge")
      result$ridge <- res.rdg[1]
    }
    if (!inherits(res.rdg, "try-error")) {
      for (iL in 1:length(ridge.lambda)) {
        if (length(ridge.lambda)>1)
          chemin<-order(abs(res.rdg$coef[,iL]),decreasing=TRUE)[1:dmax.rdg]
        if (length(ridge.lambda)==1)
          chemin<-order(abs(res.rdg$coef),decreasing=TRUE)[1:dmax.rdg]
        f.projiL <- matrix(0,nrow=n,ncol=dmax.rdg)
        r.projiL <- rep(0,dmax.rdg)
        m.rdgiL <- list(NULL)
        for (h in 1:dmax.rdg) {
          m.rdgiL[[h]] <- chemin[1:h]
          Pr <- ProjY(Y,cbind(rep(1,n),X[,chemin[1:h]]),h+1)
          f.projiL[,h] <- Pr$Proj
          r.projiL[h] <- Pr$rg
        }
        f.proj <- cbind(f.proj,f.projiL)
        r.proj <- c(r.proj,r.projiL)
        m.rdg[((iL-1)*dmax.rdg+1):(iL*dmax.rdg)]  <- m.rdgiL
      }
      f.proj <- cbind(rep(mean(Y),n),f.proj)
      r.proj <- c(1,r.proj)
      m.rdg <- c("Intercept",m.rdg)
      res1 <-  EstSelect5(Y, X, f.proj, m.rdg, r.proj, pen)
      if (long.output) {
        result$chatty$ridge<-list(support=m.rdg,crit=res1$crit)
      }
      result$ridge<-list(support=res1$mHat,crit=res1$crit[res1$lHat],fitted=res1$fHat)
      crit <- c(crit,res1$crit[res1$lHat])
      ll=ll+1
      supp[[ll]]<-res1$mHat
      meth <- c(meth,"ridge")
    }
  }
  # Random Forest procedure
  if (is.element("rF",method)) {
    if (verbose==TRUE) {
      print(paste("RANDOMFOREST PROCEDURE with option dmax=",min(rF.dmax,dmax),
                    "rF.lmtry=",rF.lmtry,collapse=" "))
    }
#    library(randomForest)
    ##note<< When method is \code{rF}, library \code{randomForest} is loaded.
    dmax.rF <- min(rF.dmax,dmax)
    m.rF <- list(NULL)
    f.proj <-NULL
    r.proj <- NULL
    nberr=0
    for (iL in 1:length(rF.lmtry)) {
      res.rF <-
        try(randomForest(X,Y,importance=TRUE,mtry=floor(p/rF.lmtry[iL])))
      if (inherits(res.rF, "try-error")) {
        nberr=nberr+1
        print(paste("error with method=rF, when calling randomForest with parameter
lmtry =",rF.lmtry[iL]))
        next
      }
      for (j in 1:2) {
        chemin <- order(res.rF$importance[,j],decreasing=TRUE)[1:dmax.rF]
        f.projiL <- matrix(0,nrow=n,ncol=dmax.rF)
        r.projiL <- rep(0,dmax.rF)
        m.rFiL <- list(NULL)
        for (h in 1:dmax.rF) {
          m.rFiL[[h]] <- chemin[1:h]
          Pr <-  ProjY(Y,cbind(rep(1,n),X[,chemin[1:h]]),h+1)
          f.projiL[,h] <- Pr$Proj
          r.projiL[h] <- Pr$rg
        }
        f.proj <- cbind(f.proj,f.projiL)
        r.proj <- c(r.proj,r.projiL)
        m.rF[((2*(iL-1)+j-1)*dmax+1):((2*(iL-1)+j)*dmax.rF)]  <- m.rFiL
      }
    }
    if (nberr==length(rF.lmtry)) {
      result$rF <- "error in randomForest for all values of ltry"
    }
    if (nberr!=length(rF.lmtry)) {
      f.proj <- cbind(rep(mean(Y),n),f.proj)
      r.proj <- c(1,r.proj)
      m.rF <- c("Intercept",m.rF)
      res1 <-  EstSelect5(Y, X, f.proj, m.rF, r.proj, pen)
      if (long.output) {
        result$chatty$rF<-list(support=m.rF,crit=res1$crit)
      }
      result$rF<-list(support=res1$mHat,crit=res1$crit[res1$lHat],fitted=res1$fHat)
      crit <- c(crit,res1$crit[res1$lHat])
      ll=ll+1
      supp[[ll]]<-res1$mHat
      meth <- c(meth, "rF")
    }
  }
  # pls procedure
  if (is.element("pls",method)) {
    if (verbose==TRUE) {
      print(paste(c("PLS PROCEDURE with options dmax=",min(pls.dmax,dmax),
           "pls.ncomp=",1:pls.ncomp),collapse=" "))
    }
#    library(pls)
    ##note<< When method is \code{pls}, library \code{pls} is loaded.
    ##details<< When method is \code{pls} or \code{ALpls}, the
    ##\code{LINselect} procedure is carried out considering the number
    ##of components in the \code{pls} method as the tuning
    ##parameter. \cr This tuning parameter varies from 1 to \code{pls.ncomp}.
    dmax.pls <- min(pls.dmax,dmax)
    m.pls <- list(NULL)
    f.proj<-NULL
    r.proj <- NULL
    res.pls <- try(plsr(Y~X,ncomp=pls.ncomp,method="simpls",validation="none"))
    if (inherits(res.pls, "try-error")) {
      print("error with method=pls, when calling plsr")
      result$pls <- res.pls[1]
    }
    if (!inherits(res.pls, "try-error")) {
      for (iL in 1:pls.ncomp) {
        chemin <- order(abs(res.pls$coef[,1,iL]),decreasing=TRUE)[1:dmax.pls]
        f.projiL <- matrix(0,nrow=n,ncol=dmax.pls)
        r.projiL <- rep(0,dmax.pls)
        m.plsiL <- list(NULL)
        for (h in 1:dmax.pls) {
          m.plsiL[[h]] <- chemin[1:h]
          Pr <-  ProjY(Y,cbind(rep(1,n),X[,chemin[1:h]]),h+1)
          f.projiL[,h] <-Pr$Proj
          r.projiL[h] <- Pr$rg
        }
        f.proj <- cbind(f.proj,f.projiL)
        r.proj <- c(r.proj, r.projiL)
        m.pls[((iL-1)*dmax.pls+1):(iL*dmax.pls)]  <- m.plsiL
      }
      f.proj <- cbind(rep(mean(Y),n),f.proj)
      r.proj <- c(1,r.proj)
      m.pls <- c("Intercept",m.pls)
      res1 <-  EstSelect5(Y, X, f.proj, m.pls, r.proj, pen)
      if (long.output) {
        result$chatty$pls<-list(support=m.pls,crit=res1$crit)
      }
      result$pls<-list(support=res1$mHat,crit=res1$crit[res1$lHat],fitted=res1$fHat)
      crit <- c(crit,res1$crit[res1$lHat])
      ll=ll+1
      supp[[ll]]<-res1$mHat
      meth <- c(meth, "pls")
    }
  }
  # Adaptive Lasso Ridge procedure
  if (is.element("ALridge",method)) {
    if (verbose==TRUE) {
      print("ADAPTIVE LASSO with WEIGHTS based on the RIDGE PROCEDURE")
      print(paste(c("with options max.steps=",max.steps,
                  "dmax=",min(ALridge.dmax, dmax),
                  "ALridge.lambda=",ALridge.lambda),collapse=" "))
    }
#    library(elasticnet)
#    library(MASS)
    ##note<< When method is \code{ALridge}, libraries \code{MASS} and \code{elasticnet} are loaded.
    m.lasso <- list(NULL)
    Nm <- 0
    f.proj<-NULL
    r.proj <- NULL
    res.rdg <- lm.ridge(Y~X,lambda=ALridge.lambda)
    if (inherits(res.rdg, "try-error")) {
      print("error with method=ALridge, when calling lm.ridge")
      result$ALridge <- res.rdg[1]
    }
    if (!inherits(res.rdg, "try-error")) {
      nberr <- 0
      for (iL in 1:length(ALridge.lambda)) {
        if (length(ALridge.lambda)>1) U <-  t(t(X)*abs(res.rdg$coef[,iL]))
        if (length(ALridge.lambda)==1) U <-  t(t(X)*abs(res.rdg$coef))
        res.lasso <-
          try(enet(U,Y,lambda=0,intercept=TRUE,normalize=FALSE,max.steps=n))
        if (inherits(res.lasso, "try-error")) {
          nberr=nberr+1
          print(paste("error with method=ALridge, when calling enet with parameter
lridge =",ALridge.lambda[iL]))
          next
        }
        BetaHat <- array(0,c(dim(res.lasso$beta.pure)[1],p))
        BetaHat[,res.lasso$allset] <- res.lasso$beta.pure
        DimBetaHat <- apply(BetaHat!=0,1,sum)
        Nmod.lasso  <- dim(BetaHat)[1]
        ff <- calc.fhat(Y,X,BetaHat,Nmod.lasso,n,p)
        f.projiL <- ff$fhat
        r.projiL <- ff$r
        m.lassoiL <- ff$m
        Nmod <-  which(DimBetaHat<=min(ALridge.dmax, dmax))
        f.proj <- cbind(f.proj,f.projiL[,Nmod])
        r.proj <- c(r.proj, r.projiL[Nmod])
        m.lasso[(Nm+1):(Nm+length(Nmod))]  <- m.lassoiL[Nmod]
        Nm<- Nm+length(Nmod)
      }
      if (nberr==length(ALridge.lambda)) {
        result$ALridge <- "error in enet for all values of lridge"
      }
      if (nberr!=length(ALridge.lambda)) {
        res1 <-  EstSelect5(Y, X, f.proj, m.lasso, r.proj, pen)
        if (long.output) {
          result$chatty$ALridge<-list(support=m.lasso,crit=res1$crit)
        }
        result$ALridge<-list(support=res1$mHat,crit=res1$crit[res1$lHat],fitted=res1$fHat)
        crit <- c(crit,res1$crit[res1$lHat])
        ll=ll+1
        supp[[ll]]<-res1$mHat
        meth <- c(meth,"ALridge")
      }
    }
  }
  # Adaptive Lasso pls procedure
  if (is.element("ALpls",method)) {
    if (verbose==TRUE) {
      print("ADAPTIVE LASSO  with WEIGHTS based on the PLS PROCEDURE")
      print(paste(c("with options max.steps=",max.steps,
                  "dmax=",min(ALpls.dmax, dmax),
                  "pls.ncomp=",1:ALpls.ncomp),collapse=" "))
     }
#    library(elasticnet)
#    library(pls)
    ##note<< When method is \code{ALpls}, libraries \code{pls} and \code{elasticnet} are loaded.
    m.lasso <- list(NULL)
    Nm <- 0
    f.proj<-NULL
    r.proj <- NULL
    res.pls <-
      try(plsr(Y~X,ncomp=ALpls.ncomp,method="simpls",validation="none"))
    if (inherits(res.pls, "try-error")) {
      print("error with method=ALpls, when calling plsr")
      result$ALpls <- res.pls[1]
    }
    if (!inherits(res.pls, "try-error")) {
      nberr <- 0
      for (iL in 1:ALpls.ncomp) {
        U <-  t(t(X)*abs(res.pls$coef[,1,iL]))
        res.lasso <-
          try(enet(U,Y,lambda=0,intercept=TRUE,normalize=FALSE,max.steps=n))
        if (inherits(res.lasso, "try-error")) {
          nberr=nberr+1
          print(paste("error with method=ALpls, when calling enet with",
                      iL,"pls components"))
          next
        }
        if (!inherits(res.lasso, "try-error")) {
          BetaHat <- array(0,c(dim(res.lasso$beta.pure)[1],p))
          BetaHat[,res.lasso$allset] <- res.lasso$beta.pure
          DimBetaHat <- apply(BetaHat!=0,1,sum)
          Nmod.lasso  <- dim(BetaHat)[1]
          ff <- calc.fhat(Y,X,BetaHat,Nmod.lasso,n,p)
          f.projiL <- ff$fhat
          r.projiL <- ff$r
          m.lassoiL <- ff$m
          Nmod <-  which(DimBetaHat<=min(ALpls.dmax, dmax))
          f.proj <- cbind(f.proj,f.projiL[,Nmod])
          r.proj <-  c(r.proj,r.projiL[Nmod])
          m.lasso[(Nm+1):(Nm+length(Nmod))]  <- m.lassoiL[Nmod]
          Nm<- Nm+length(Nmod)
        }
      }
      if (nberr==ALpls.ncomp) {
        result$ALpls<-"error in enet for all values of the number of components"
      }
      if (nberr!=length(ALpls.ncomp)) {
        res1 <-  EstSelect5(Y, X, f.proj, m.lasso, r.proj, pen)
        if (long.output) {
          result$chatty$ALpls<-list(support=m.lasso,crit=res1$crit)
        }
        resultALpls<-list(support=res1$mHat,crit=res1$crit[res1$lHat],fitted=res1$fHat)
        crit <- c(crit,res1$crit[res1$lHat])
        ll=ll+1
        supp[[ll]]<-res1$mHat
        meth <- c(meth, "ALpls")
      }
    }
  }
  # Exhaustive procedure
  if (is.element("exhaustive",method)) {
    if (verbose==TRUE) print("EXHAUSTIVE PROCEDURE")
#    library(gtools)
    ##note<< When method is \code{exhaustive}, library \code{gtools}
    ##is loaded.
    ##details<< When method is \code{exhaustive}, the maximum
    ##number of variate d is calculated as
    ##follows.\cr
    ##Let q be the largest integer such that \code{choose(p,q)} <
    ##\code{exhaustive.maxdim}. Then d = \code{min(q, exhaustive.dmax,dmax)}.
    dmax.ex1 <- min(exhaustive.dmax,dmax)
    dim <- which(choose(p,1:p) > exhaustive.maxdim)
    if (length(dim)==0) dmax.ex <- min(dmax.ex1,p)
    if (length(dim)>0) dmax.ex <- min(dim[1]-1,dmax.ex1)
    if (verbose==TRUE) {
      print(paste("dmax, the maximum number of selected variates equals",dmax.ex))
      if ((length(dim)>0)&(dmax.ex<dmax.ex1)) {
        print("to increase dmax, choose a larger value for the parameter exhaustive.maxdim")
      }
    }
    m.ex <- list(NULL)
    f.proj <- NULL
    r.proj <- NULL
    Nm <- 0
# SH : 21/02/2013 boucle a ecrire en C
    for (h in 1:dmax.ex) {
      nmod <- choose(p,h)
      if (verbose==TRUE) {
        print(paste("selection of the best",h,"variates"))
        print(paste("the number of subsets of",h,
                    "variates among",p,"equals",nmod))
      }
      mod <- combinations(p,h)
      f.projh <- matrix(0,nrow=n,ncol=nmod)
      r.projh <- rep(0,nmod)
      m.exh <- list(NULL)
      for (iL in (1:nmod)) {
        Pr <-  ProjY(Y,cbind(rep(1,n),X[,mod[iL,]]),h+1)
        f.projh[,iL] <- Pr$Proj
        r.projh[iL] <- Pr$rg
        m.exh[[iL]] <- mod[iL,]
      }
      f.proj <- cbind(f.proj,f.projh)
      r.proj <- c(r.proj,r.projh)
      m.ex[(Nm+1):(Nm+length(m.exh))] <- m.exh
      Nm <- Nm+length(m.exh)
    }
    f.proj <- cbind(rep(mean(Y),n),f.proj)
    r.proj <- c(1,r.proj)
    m.ex <- c("Intercept",m.ex)
    res1 <- EstSelect5(Y, X, f.proj, m.ex, r.proj, pen)
    if (long.output) {
      result$chatty$exhaustive<-list(support=m.ex,crit=res1$crit)
    }
    result$exhaustive<-list(support=res1$mHat,crit=res1$crit[res1$lHat],fitted=res1$fHat)
    crit <- c(crit,res1$crit[res1$lHat])
    ll=ll+1
    supp[[ll]]<-res1$mHat
    meth <- c(meth, "exhaustive")
  }
  ii <- which.min(crit)
  result$summary <- list(support=supp[ii],
                         crit=crit[ii],
                         method=meth[which(crit[ii]==crit)])
                             
  return(result)
  ### A list with at least \code{length(method)}
  ### components. \cr 
  ### For each procedure in \code{method}  a list with components \cr
  ### \itemize{
  ###    \item{\code{support}: vector of integers. Estimated support of the
  ###              parameters \eqn{\beta} for the considered procedure.}
  ###   \item{\code{crit}: scalar equals to the LINselect criteria
  ###   calculated in the estimated support.}
  ###   \item{\code{fitted}: vector  with length n. Fitted value of
  ###       the response calculated when the support of \eqn{ \beta}
  ###       equals \code{support}.}
  ###   \item{\code{coef}: vector whose first component is the estimated
  ###             intercept. \cr The other components are the estimated non zero
  ###             coefficients when the support of \eqn{ \beta}
  ###       equals \code{support}.}
  ### }
  ### If \code{length(method)} > 1, the additional component \code{summary} is a list with three
  ### components:
  ### \itemize{ 
  ###    \item{\code{support}: vector of integers. Estimated support of the
  ###              parameters \eqn{\beta} corresponding to the minimum
  ###              of the criteria among all procedures.}
  ###    \item{\code{crit}:  scalar. Minimum value of the
  ###    criteria among all procedures.} 
  ###     \item{\code{method}: vector of characters. Names of the
  ###     procedures  for
  ###         which the minimum is reached}
  ### }
  ### If \code{pen.crit = NULL}, the component \code{pen.crit} gives the
  ### values of the penalty calculated by the function \code{penalty}.
  ### If \code{long.output} is TRUE the component named
  ### \code{chatty} is a list  with \code{length(method)}
  ### components. \cr 
  ### For each procedure in \code{method}, a list with components
  ### \itemize{
  ###   \item{\code{support} where \code{support[[l]]} is a vector of
  ###   integers containing an estimator of the support of the
  ###   parameters \eqn{ \beta}.}
  ###   \item{\code{crit} : vector where \code{crit[l]} contains the
  ###   value of the LINselect criteria calculated in
  ###   \code{support[[l]]}.}
  ### }

}#  fin VARselect
