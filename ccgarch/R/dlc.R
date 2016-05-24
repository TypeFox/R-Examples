# the gradient and Hessian of the DCC part of the log-likelihood function
dlc <- function(dcc.para, B, u, h, model){
   nobs <- dim(h)[1]
   ndim <- dim(h)[2]
   dccpar1 <- dcc.para[1]
   dccpar2 <- dcc.para[2]
   hlag <- rbind(colMeans(h), h)
   dhw <- vec.garch.derivative(u, B, h)
   dhwlag <- rbind(colMeans(dhw), dhw)
      inv.sq.h <- 1/sqrt(h)
      inv.sq.h.lag <- rbind(colMeans(inv.sq.h), inv.sq.h)
      z <- u*inv.sq.h
      zlag <- rbind(colMeans(z), z)
   Q <- cov(z)
   dcc <- dcc.est(z, dcc.para)
   P <- dcc$DCC
   Plag <- rbind(colMeans(P), P)
   Qt <- dcc$Q
   Qtlag <- rbind(colMeans(Qt), Qt)
   In <- diag(ndim)
      ind <- as.vector(rbind(1, In, In))
      if(model=="diagonal"){
         npar.h <- 3*ndim
      } else {
         npar.h <- ndim*(2*ndim+1)
      }
   Qtilde <- 1/sqrt(Qt[,as.vector(In)==1]) 
   dvecQ <- matrix(0, nobs, 2*ndim^2)
   dvecP <- matrix(0, nobs, 2*ndim^2)
   dlc <- matrix(0, nobs, 2)
   d2lc <- matrix(0, nobs, 4)
   const <- -rbind(as.vector(Q), as.vector(Q))
   dvecQt <- matrix(0, 2, ndim^2)
   dwdvecQt <- matrix(0, npar.h, ndim^2)
   dwdvecQ <- matrix(0, nobs, npar.h*ndim^2)
   dfdwd2lc <- matrix(0, nobs, 2*npar.h)
   for ( i in 1:nobs){
      dvecQt <- const + rbind(as.vector(outer(zlag[i,], zlag[i,])), Qtlag[i,]) + dccpar2*dvecQt
      dvecQ[i,] <- as.vector(dvecQt)
      invQtilde <- diag(1/Qtilde[i,])
      Pt <- matrix(P[i,], ndim, ndim)
      invPt <- solve(Pt)
      Z <- Pt%x%invQtilde; tZ <- invQtilde%x%Pt
      QPQ <- invQtilde%x%invQtilde - 0.5*diag(as.vector(invQtilde))%*%( Z + tZ )
      dvecPt <- dvecQt%*%QPQ
      dvecPt[,as.vector(In)==1] <- 0
      dvecP[i,] <- as.vector(dvecPt)
      zz <- outer(z[i,], z[i,])
      dlc[i,] <- -0.5*as.vector(dvecPt%*%as.vector(invPt - invPt%*%zz%*%invPt))
   
      d2lc[i,] <- 0.5*as.vector(dvecPt%*%(invPt%x%invPt)%*%t(dvecPt))
         dhwtlag <- matrix(dhwlag[i,], ncol=ndim)
         dhwt <- matrix(dhw[i, ], ncol=ndim)
         if(model=="diagonal"){
            dhwtlag <- dhwtlag[ind==1,]
            dhwt <- dhwt[ind == 1, ]
         }
         dwdvecVtlag <- matrix(0, npar.h, ndim^2)
         dwdvecVt <- matrix(0, npar.h, ndim^2)
         dwdvecVtlag[,as.vector(In)==1] <- dhwtlag
         dwdvecVt[, as.vector(In) == 1] <- dhwt
         Dlag <- diag(hlag[i,])
         Dt <- diag(sqrt(h[i,]))
         tmpDIn <- diag(Dt%x%In + In%x%Dt)
         tmpDInlag <- diag(Dlag%x%In + In%x%Dlag)
         dwdvecDlag <- dwdvecVtlag%*%diag(1/tmpDInlag)
         dwdvecD <- dwdvecVtlag%*%diag(1/tmpDIn)
         invD <- diag(inv.sq.h[i,])
         Ptlag <- matrix(Plag[i,], ndim, ndim)
         invDlag <- diag(inv.sq.h.lag[i,])
         
         ZD <- Ptlag%x%invDlag
         DZ <- invDlag%x%Ptlag
         dwdvecQt <- -0.5*dccpar1*dwdvecDlag%*%(ZD + DZ) + dccpar2*dwdvecQt
         dwdvecQ[i,] <- as.vector(dwdvecQt)
         dwdvecPt <- t(dwdvecQt%*%QPQ)
         
         dfdwd2lc[i,] <- as.vector(0.5*dvecPt%*%(invPt%x%invPt)%*%dwdvecPt + 0.5*dvecPt%*%((invPt%*%invD)%x%In + In%x%(invPt%*%invD))%*%t(dwdvecD))
   }
   list(dlc=dlc, dvecP = dvecP, dvecQ = dvecQ, d2lc=d2lc, dfdwd2lc=dfdwd2lc)
}

