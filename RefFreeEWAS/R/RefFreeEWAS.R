#################################################################################
# RefFreeEwasModel: Reference-free cell-mixture-adjusted EWAS
#################################################################################
RefFreeEwasModel <- function(
  Y,   
  X,   
  K,   
  smallOutput=FALSE  
){
  n1 <- dim(Y)[1]
  n2 <- dim(X)[1]
  pdim <- dim(X)[2]

  noMiss <- !apply(is.na(Y),1,any)
  allObs <- all(noMiss)

  HX <- solve(t(X)%*%X)
  PX <- X %*% HX

  if(allObs){ # If no missing value
    Bstar <- Y %*% PX
 
    muStar <- Bstar %*% t(X)
    Estar <- Y - muStar

    sstar <- sqrt(apply(Estar*Estar,1,sum)/(n2-pdim))

    BetaE <- cbind(Bstar, Estar)
    svdStar <- svdSafe(BetaE)

    Lambda <- t(svdStar$d[1:K]*t(svdStar$u[,1:K]))
    U <- svdStar$v[,1:K]
  }
  else{ #If missing values, do as much as possible on one fell swoop,
        # and for the rest, do on a CpG-by-CpG basis
    Bstar <- matrix(NA, n1, pdim)
    degfree <- rep(n2-pdim, n1)
    Bstar[noMiss,] <- Y[noMiss,] %*% PX
    whichMiss <- which(!noMiss)
    nMiss <- length(whichMiss)

    for(j in 1:nMiss){
       jj <- whichMiss[j]
       mflag <- !is.na(Y[jj,])
       Xj <- X[mflag,,drop=FALSE]
       HXj <- solve(t(Xj)%*%Xj)
       PXj <- Xj %*% HXj
       Bstar[jj,] <- Y[jj,mflag,drop=FALSE]%*%PXj
       degfree[jj] <- sum(mflag)-pdim
    }

    muStar <- Bstar %*% t(X)
    Estar <- Y - muStar

    sstar <- sqrt(apply(Estar*Estar,1,sum,na.rm=TRUE)/degfree)
    BetaE <- cbind(Bstar, Estar)
    svdStar <- svdSafe(BetaE[noMiss,])

    Lambda <- matrix(NA, n1, K)
    Lambda[noMiss,] <- t(svdStar$d[1:K]*t(svdStar$u[,1:K]))
    U <- svdStar$v[,1:K]
    for(j in 1:nMiss){
       jj <- whichMiss[j]
       mflag <- c(rep(TRUE,pdim), !is.na(Y[jj,]))
       Uj <- U[mflag,,drop=FALSE]
       Lambda[jj,] <- solve(t(Uj)%*%Uj,t(Uj)%*%BetaE[jj,mflag])
    }
  }

  LambdaProjBstar <- solve(t(Lambda)%*%Lambda, t(Lambda)%*%Bstar)
  Beta <- Bstar - Lambda %*% LambdaProjBstar

  out <- list(Bstar=Bstar, Beta=Beta, sstar=sstar, Lambda=Lambda, U=U, d=svdStar$d)

  if(!smallOutput) {
     muStar <- ifelse(muStar<0.00001,0.00001,muStar)
     muStar <- ifelse(muStar>0.99999,0.99999,muStar)
     out$dispersion <- sqrt(muStar*(1-muStar))
     out$E <- Estar/out$dispersion
     out$X <- X
  }
  out
  class(out) <- "RefFreeEwasModel"
  out
}

print.RefFreeEwasModel <- function(x,...){
   cat("Reference Free EWAS Model\n\n")
   cat("Assay matrix: ", dim(x$Beta)[1], " features\n")
   if(!is.null(x$X)) cat("Design matrix: ", dim(x$X)[1], 
      " subjects x ", dim(x$X)[2]," covariates\n\n", sep="")
   else cat("(small output version)\n\n")
 
   cat(dim(x$Lambda)[2]," latent variables\n\n")
}


#################################################################################
# BootOneRefFreeEwasModel: Create *one( bootstrap sample from
#     reference-free cell-mixture-adjusted EWAS
#################################################################################
BootOneRefFreeEwasModel <- function(mod){
  n2 <- dim(mod$X)[1]
  iboot <- sample(1:n2, n2, replace=TRUE)

  mu <- mod$Bstar %*% t(mod$X)

  return(mu + mod$dispersion*mod$E[,iboot])
}



#################################################################################
# BootRefFreeEwasModel: Create several bootstrap samples from
#     reference-free cell-mixture-adjusted EWAS
#     and return Beta estimates
#################################################################################
BootRefFreeEwasModel <- function(
   mod,   
   nboot  
){
   BetaBoot <- array(NA, dim=c(dim(mod$Beta),2,nboot))
   dimnames(BetaBoot)[1:2] <- dimnames(mod$Beta) 
   dimnames(BetaBoot)[[3]] <- c("B","B*")
   dimnames(BetaBoot)[[4]] <- 1:nboot
   attr(BetaBoot,"nSample") <- dim(mod$X)[1]

   for(r in 1:nboot){
      isError <- TRUE
      while(isError){
        catchError <- try({
          Yboot <- BootOneRefFreeEwasModel(mod)
          bootFit <- RefFreeEwasModel(Yboot, mod$X, 
            dim(mod$Lambda)[2], smallOutput=TRUE)
          BetaBoot[,,1,r] <- bootFit$Beta
          BetaBoot[,,2,r] <- bootFit$Bstar
        })
        isError <- inherits(catchError,"try-error")
      }
      if(r%%10==0) cat(r,"\n")
   }
   class(BetaBoot) <- "BootRefFreeEwasModel"
   BetaBoot
}

summary.BootRefFreeEwasModel <- function(object,...){
  
   x <- object
 
   out <- array(NA, dim=c(dim(x)[1:3],2))
   dimnames(out)[1:3] <- dimnames(x)[1:3]
   dimnames(out)[[4]] <- c("mean","sd")

   out[,,,1] <- apply(x,c(1:3),mean)
   out[,,,2] <- apply(x,c(1:3),sd)

   class(out) <- "summaryBootRefFreeEwasModel"
   attr(out,"nBoot") <- dim(x)[4]
   attr(out,"nSample") <- attr(x,"nSample")

   out
}

print.summaryBootRefFreeEwasModel <- function(x,...){
   cat(attr(x,"nBoot"),"bootstrap samples, n =", attr(x,"nSample"),"subjects\n\nBeta Mean\n")
   print(x[1:6,,1,1],...)
   cat("\nBeta Standard Deviation\n")
   print(x[1:6,,1,2],...)
}

print.BootRefFreeEwasModel <- function(x,...){
   print(summary(x),...)
}

#################################################################################
# svdSafe: svd that traps errors and switches to QR when necessary
#################################################################################
svdSafe <- function(X){
  sv <- try(svd(X), silent=TRUE)
  if(inherits(sv,"try-error")){
    warning("SVD algorithm failed, using QR-decomposition instead")
    QR <- qr(X)
    sv <- list(d=rep(1,dim(X)[2]))
    sv$u <- qr.Q(QR)
    sv$v <- t(qr.R(QR))
  }
  sv
}

#################################################################################
# EstDimIC: Estimate dimension using AIC and BIC
#################################################################################
EstDimIC <- function (
  Rmat,             
  Krange=0:25
){
  N1 <- dim(Rmat)[1]
  N2 <- dim(Rmat)[2]
  svdRmat <- svdSafe(Rmat)
  nK = length(Krange)
  tmpAIC <- tmpBIC <- rep(NA,nK)
  for(Ktest in Krange){
    if(Ktest==0) tmpRminLU <- Rmat
    else tmpRminLU <- (Rmat -
        svdRmat$u[,1:Ktest] %*% (svdRmat$d[1:Ktest] * t(svdRmat$v[,1:Ktest])))
    tmpSigSq <- apply(tmpRminLU*tmpRminLU,1,sum)/N2
    tmpAIC[Ktest+1] <- 2*(N1+Ktest*(N1+N2)) + N1*N2 + N2*sum(log(tmpSigSq))
    tmpBIC[Ktest+1] <- log(N2)*(N1+Ktest*(N1+N2)) + N1*N2 + N2*sum(log(tmpSigSq))
  }
  list(icTable=cbind(K=Krange,AIC=tmpAIC,BIC=tmpBIC), 
       best=Krange[c(AIC=which.min(tmpAIC),BIC=which.min(tmpBIC))])
}

#################################################################################
# omnibusBoot: Omnibus test of association using bootstraps
#################################################################################
omnibusBoot <- function(est, boots, denDegFree){
  nFeature <- length(est)
  se <- apply(boots,1,sd)
  pv <- 2*pt(-abs(est)/se,denDegFree)
  pvNull <- 2*pt(-(1/se)*abs(sweep(boots,1,apply(boots,1,mean),"-")),denDegFree+1)

  ks <- max(abs(sort(pv)-(1:nFeature-0.5)/nFeature))
  ksNull <-apply(pvNull,2,function(u)max(abs(sort(u)-(1:nFeature-0.5)/nFeature)))
  mean(ks<=ksNull)
}

