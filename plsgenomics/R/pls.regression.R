### pls.regression.R  (2005-05-10)
###
###     Multivariate Partial Least Squares Regression
###
### Copyright 2004-05 Anne-Laure Boulesteix and Korbinian Strimmer
###
### Part of the code was adopted from the pls.pcr package by Ron Wehrens
###
###
### This file is part of the `plsgenomics' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA





#
# original simpls function from "pls.pcr" package
# modified so that weights R are also returned
# 

standard.simpls <- function (Xtrain, Ytrain, Xtest=NULL, ncomp=NULL)
{
    X<-scale(Xtrain,center=TRUE,scale=FALSE)
    meanX<-attributes(X)$"scaled:center"
    if (is.vector(Ytrain))
     {
     Ytrain<-matrix(Ytrain,length(Ytrain),1)
     }
    Y<-scale(Ytrain,center=TRUE,scale=FALSE)
    n <- dim(X)[1]
    p <- dim(X)[2]
    m <- dim(Y)[2]
    if (is.null(ncomp))
     {
     ncomp<-min(n,p)
     }
     
    S <- crossprod(X, Y)
    RR <- matrix(0, ncol = max(ncomp), nrow = p)
    PP <- matrix(0, ncol = max(ncomp), nrow = p)
    QQ <- matrix(0, ncol = max(ncomp), nrow = m)
    TT <- matrix(0, ncol = max(ncomp), nrow = n)
    VV <- matrix(0, ncol = max(ncomp), nrow = p)
    UU <- matrix(0, ncol = max(ncomp), nrow = n)
    B <- array(0, c(dim(X)[2], dim(Y)[2], length(ncomp)))
    if (!is.null(Xtest))
        {
	if (is.vector(Xtest))
	 {
	 Xtest<-matrix(Xtest,1,length(Xtest))
	 }
        Ypred <- array(0, c(dim(Xtest)[1], m, length(ncomp)))
	}
    for (a in 1:max(ncomp)) {
        qq <- svd(S)$v[, 1]
        rr <- S %*% qq
        tt <- scale(X %*% rr, scale = FALSE)
        tnorm <- sqrt(sum(tt * tt))
        tt <- tt/tnorm
        rr <- rr/tnorm
        pp <- crossprod(X, tt)
        qq <- crossprod(Y, tt)
        uu <- Y %*% qq
        vv <- pp
        if (a > 1) {
            vv <- vv - VV %*% crossprod(VV, pp)
            uu <- uu - TT %*% crossprod(TT, uu)
        }
        vv <- vv/sqrt(sum(vv * vv))
        S <- S - vv %*% crossprod(vv, S)
        RR[, a] <- rr
        TT[, a] <- tt
        PP[, a] <- pp
        QQ[, a] <- qq
        VV[, a] <- vv
        UU[, a] <- uu
        if (!is.na(i <- match(a, ncomp))) {
            B[, , i] <- RR[, 1:a, drop = FALSE] %*% t(QQ[, 1:a,
                drop = FALSE])
            if (!is.null(Xtest))
	     {
	     Xtest<-scale(Xtest,scale=FALSE,center=meanX)
             Ypred[, , i] <- Xtest %*% B[, , i]
	     }
        }
    }
    if (length(ncomp)==1)
     {
     B<-B[,,1]
     }
    if (!is.null(Xtest))
        list(B = B, Ypred = Ypred, P = PP, Q = QQ, T = TT, R=RR, meanX=meanX)
    else list(B = B, P = PP, Q = QQ, T = TT, R=RR, meanX=meanX)
}


# 
# original simpls function returns orthonormal Xscores                                      
# but the weight vectors r_i are NOT standardized to length 1.
#
# instead, this functions returns orthogonal Xscores
# with weight vectors r_i of length 1
#                             

unitr.simpls <- function (Xtrain, Ytrain, Xtest=NULL, ncomp=NULL)
{
  pls.out <- standard.simpls(Xtrain=Xtrain, Ytrain=Ytrain, Xtest=Xtest, ncomp)

  # norm of vector
  euclidian.norm <- function(xvec)
  {
    return( sqrt(sum(xvec*xvec)) )
  }

  # Compute norm of weight vectors r_i
  R.norm <- apply(pls.out$R, 2, euclidian.norm)
  
 
  # Scaling matrices
  if (length(R.norm)==1)
   {
   M<-matrix(R.norm,1,1)
   Mi<-matrix(1/R.norm,1,1)
   }
  else
   {  
   M <- diag(R.norm)
   Mi <- diag(1/R.norm)
   }
  # Transform output matrices
  Rnew <- pls.out$R %*% Mi
  Tnew <- pls.out$T %*% Mi
  
  Qnew <- pls.out$Q %*% M
  Pnew <- pls.out$P %*% M

  # B and Ypred are invariant !


  if (!is.null(Xtest))
        list(B = pls.out$B, Ypred = pls.out$Ypred, P = Pnew, Q = Qnew,
             T = Tnew, R=Rnew, meanX=pls.out$meanX)
    else list(B = pls.out$B, P = Pnew, Q = Qnew, T = Tnew, R=Rnew,
    meanX=pls.out$meanX)
}

######################

pls.regression <- function(Xtrain, Ytrain, Xtest=NULL, ncomp=NULL,  unit.weights=TRUE)
{
  if (unit.weights==TRUE)
  {
    return( unitr.simpls(Xtrain, Ytrain, Xtest=Xtest, ncomp) )
  }
  else
  {
    return( standard.simpls(Xtrain, Ytrain, Xtest=Xtest, ncomp) )
  }
}

