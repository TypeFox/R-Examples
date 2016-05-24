coxphw.fit <- function
(
 obj,
 id,
 weights,
 PARMS,
 CARDS=NULL ,
 sorted=FALSE,
 pc=TRUE,
 pc.time=TRUE,
 standardize=TRUE,
 fixed=NULL
 )
### fitter function
### 2010-07-07
{       robust <- (PARMS[3]==1 | PARMS[3]==3)
        jack <- (PARMS[3]>1)
        k <- ncol(obj$mm1)    # number covariates w/o time-dep effects
        k2 <- k + obj$NTDE
        NTDE <- obj$NTDE
        if(NTDE>0 & k>1) pc<-FALSE
        kk <- k
        maxid <- max(id)
        
        ## standardize model matrix (next 5 lines old code 2010-07)
#        sd1 <- sd(obj$mm1)
#        sd2 <- sd(obj$timedata)
#        Z.sd <- c(sd1, sd2 * sd1[obj$timeind])
#        ZxZ <- as.matrix(Z.sd) %*% t(as.matrix(Z.sd)) # used often to restandardize ...
#        obj$mm1 <- scale(obj$mm1, FALSE, sd1)
       ## next lines new code: 
        if (pc.time==TRUE & ncol(obj$timedata)>1) {
#       pc<-FALSE
            pc.t1<-prcomp(obj$timedata, center=FALSE)
            rot.mat.t<-pc.t1$rotation
            obj$timedata<-predict(pc.t1,obj$timedata)
        }
        mm1.orig<-obj$mm1
        if (pc==TRUE & sum(obj$ind[1:kk])>1) {
                   pc1<-prcomp(obj$mm1, center=FALSE)
                   rot.mat<-pc1$rotation
                   obj$mm1<-predict(pc1,obj$mm1)
               }
       if(standardize){
          sd2 <- apply(as.matrix(obj$timedata),2,sd)
          sd1 <- apply(as.matrix(obj$mm1),2,sd)
          sd1.orig <- apply(as.matrix(mm1.orig),2,sd)
       }
       else {
          sd2<-rep(1,ncol(obj$timedata))
          sd1 <- rep(1, ncol(obj$mm1))
          sd1.orig <- rep(1, ncol(mm1.orig))
       }
       Z.sd <- c(sd1, sd2 * sd1.orig[obj$timeind])
       ZxZ <- as.matrix(Z.sd) %*% t(as.matrix(Z.sd))
       obj$mm1 <- scale(obj$mm1, FALSE, sd1)
       mm1.orig <- scale(mm1.orig, FALSE, sd1.orig)
       ### end new code
       if(!is.null(fixed)) fixed[!is.na(fixed)]<-fixed[!is.na(fixed)]*Z.sd[!is.na(fixed)]
        
        ##if(ind.offset)
        obj$mm1o <- if(PARMS[16] != 0) cbind(obj$offset.values, obj$mm1) else obj$mm1
        obj$timedata <- scale(obj$timedata, FALSE, sd2)
        mmm <- cbind(obj$mm1, obj$timedata) # model matrix inc. time data
        if(is.null(CARDS))
          CARDS <- cbind(obj$mm1o, obj$resp, weights, obj$timedata, id)
        if(!sorted) CARDS <- CARDS[order(obj$resp[, 2], -obj$resp[, 3]), ]
        ##   if (offset) {
        ##    IOARRAY[1,1]<-0    # first variable is offset
        ##    IOARRAY[2,1]<-Z.sd[1]    # first variable is offset
        ##   }
        DFBETA <- matrix(0, maxid, k2)
        IOARRAY <- rbind(rep(1, k2), matrix(0, 2+3*k2, k2))       #changed georg 090604
        if(!is.null(fixed)){
          IOARRAY[1,!is.na(fixed)]<-0
          IOARRAY[2,]<-fixed
          IOARRAY[2,is.na(fixed)]<-0
        }
        if(obj$NTDE >0 )
          IOARRAY[4, (k+1):k2] <- obj$timeind
        
        ## --------------- Aufruf Fortran-Routine WEIGHTEDCOX ------------
        storage.mode(CARDS) <- storage.mode(PARMS) <- storage.mode(IOARRAY) <- storage.mode(DFBETA) <- "double"
        value <- .Fortran("weightedcox",
                          cards=CARDS,
                          outpar = PARMS,
                          outtab = IOARRAY,
                          dfbetaresid = DFBETA,
                          PACKAGE="coxphw")
        if(value$outpar[8])
          warning("Error in Fortran routine weightedox; parms8 <> 0")
  #      std.cov<-matrix(value$outtab[4:(k2+3), ], ncol=k2) / ZxZ
  #      coefs<-value$outtab[3,  ] / Z.sd
  #      wald.chisq<-t(coefs)%*%solve(std.cov)%*%coefs
  #      value$outpar[9]<-wald.chisq
#        coefs <- value$coefs / Z.sd
        coef.orig <- value$outtab[3,  ]
        coefs <- coef.orig / Z.sd
        dfbeta.resid=NULL
        cov.ls <- matrix(value$outtab[4:(k2+3), ], ncol=k2) / ZxZ # is matrix() needed?..
        if(robust) cov.lw <- matrix(value$outtab[(k2+4):(3+2*k2), ], ncol=k2) / ZxZ
        else cov.lw<- NULL
        if(jack) cov.j <- matrix(value$outtab[(2*k2+4):(3+3*k2), ], ncol=k2) / ZxZ
        else cov.j <- NULL
        if(robust | jack) dfbeta.resid <- value$dfbetaresid / matrix(Z.sd, maxid, k2, byrow=TRUE)
        if (pc==TRUE & sum(obj$ind[1:kk])>1) {
           rot.mat.blowup<-matrix(0,k+NTDE,k+NTDE)
           diag(rot.mat.blowup)<-1
           rot.mat.blowup[(1:k),(1:k)]<-rot.mat
           coefs <- rot.mat.blowup %*% coefs
           cov.ls <- rot.mat.blowup %*% cov.ls  %*%  t(rot.mat.blowup)
           if(robust) {
              cov.lw <-rot.mat.blowup %*% cov.lw  %*%  t(rot.mat.blowup)
              if(!jack) dfbeta.resid <- t(rot.mat.blowup %*% t(dfbeta.resid))
              }
           if(jack) {
             cov.j <-rot.mat.blowup %*% cov.j  %*%  t(rot.mat.blowup)
              dfbeta.resid <- t(rot.mat.blowup %*% t(dfbeta.resid))
            }
             
          }
        if (pc.time==TRUE &  ncol(obj$timedata)>1) {
           rot.mat.t.blowup<-matrix(0,k+NTDE,k+NTDE)
           diag(rot.mat.t.blowup)<-1
           rot.mat.t.blowup[(k+1):(k+NTDE),(k+1):(k+NTDE)]<-rot.mat.t
           coefs <- rot.mat.t.blowup %*% coefs
           cov.ls <- rot.mat.t.blowup %*% cov.ls %*% t(rot.mat.t.blowup)
           if(robust) {
             cov.lw <- rot.mat.t.blowup %*% cov.lw %*% t(rot.mat.t.blowup)
              if(!jack) dfbeta.resid <- t(rot.mat.t.blowup %*% t(dfbeta.resid))
              }
           if(jack) {
              cov.j <- rot.mat.t.blowup %*% cov.j %*% t(rot.mat.t.blowup)
              dfbeta.resid <- t(rot.mat.t.blowup %*% t(dfbeta.resid))
            }
              
        }

        res <- list(
                    cards=value$cards,
                    outpar=value$outpar,
                    outtab=matrix(value$outtab, nrow=3+3*k2),   # changed georg 090604
                    dfbeta.resid=dfbeta.resid,
                    coef.orig=value$outtab[3,  ], 
                    coefs=coefs, # coefficients
                    cov.ls=cov.ls,      # covariances
                    cov.lw=cov.lw,
                    cov.j=cov.j,
                    Z.sd=Z.sd,
                    ZxZ=ZxZ,
                    abs.score=value$outpar[7],
                    mmm=mmm             # model matrix
                    )
        res
}
