#******************************************************************************* 
#
# Estimation for Multivariate Normal Data with Monotone Missingness
# Copyright (C) 2007, University of Cambridge
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (bobby@statslab.cam.ac.uk)
#
#*******************************************************************************


## check.QP:
##
## check the Quadratic Programming arguments, as described
## by solve.QP in library(quadprog) as encoded in QP
## and reorder by nao

check.QP <- function(QP, m, nao, oo)
  {
    ## check if we should return the null default solve.QP params
    if(is.null(QP) || (is.logical(QP) && QP == FALSE))
      return(list(cols=NULL, dvec=NULL, dmu=FALSE, Amat=NULL, b0=NULL,
                  mu.constr=0, q=0, meq=0, o=NULL))
    
    ## use the default.QP, and then continue to check
    if(is.logical(QP) && length(QP) == 1 && QP == TRUE) QP <- default.QP(m)
    else if(is.numeric(QP)) { ## treat QP as a number of factors
      if(length(QP) != 1 || QP < 1 || QP >= m)
        stop("scalar QP should be an integer scalar 0 < QP < m")
      QP <- default.QP(m-QP)
    }
    
    ## make sure QP is a list
    if(!is.list(QP)) stop("QP should be a list object")

    ## turn nao into a non-null vector
    if(is.null(nao)) oo <- nao <- 1:m
      
    ## see how many initial cols of mu and Sigma are treated as factors
    if(m < QP$m) stop("should have QP$m <= ncol(y)=m")

    ## calculate the number of factors, 
    ## then just use QP's version of m
    numfact <- m - QP$m
    m <- QP$m

    ## calculate the indicies of the factors under nao
    if(numfact > 0) {
      facts <- oo[1:numfact]

      ## now get indices of the non-factors
      notfacts <- oo[-(1:numfact)]
      ## and then in monomvn ordering
      keep <- sort(notfacts)
            
      ## now calculating the re-ordering of the rows of Amat, etc
      nao <- rank(setdiff(nao,1:numfact))
    } else keep <- 1:m
    
    ## check the Amat argument
    Amat <- as.matrix(QP$Amat)
    if(nrow(Amat) != m) stop("QP$Amat should have", m, " rows")
    Amat <- Amat[nao,]

    ## check the b0 argument
    b0 <- as.vector(QP$b0)
    if(length(b0) != ncol(Amat))
      stop("should have length(QP$b0) == ncol(QP$Amat)")

    ## check the meq argument
    meq <- as.integer(QP$meq)
    if(length(meq) != 1 || meq > ncol(Amat) || meq < 0)
      stop("must have 0 <= QP$meq <= ncol(QP$Amat) =", ncol(Amat))
    
    ## check the dvec argument
    dvec <- as.vector(QP$dvec)
    if(is.null(dvec) || length(dvec) != m)
        stop("QP$dvec) must be a vector of length", m)

    ## check the dmu argument
    dmu <- QP$dmu
    if(is.null(dmu) || !is.logical(dmu) || length(dmu) != 1)
      stop("QP$dmu must be a scalar logical")

    ## check the mu.constr argument
    mu.constr <- QP$mu.constr
    if(is.null(mu.constr) || !is.numeric(mu.constr) ||
       any(mu.constr[-1] < 1) || length(mu.constr)-1 != mu.constr[1] ||
       any(duplicated(mu.constr[-1])) )
      stop("QP$mu.costr must be a positive integer vector with\n",
           "\tlength(mu.constr)-1 = mu.constr[1] and no duplicated\n",
           "\tentries in mu.constr[-1]")

    ## return the list
    return(list(m=m, cols=keep, dvec=dvec, dmu=dmu, Amat=Amat, b0=b0,
                mu.constr=mu.constr, q=ncol(Amat), meq=meq, o=nao))
  }


## default.QP:
##
## create the devault solve.QP setup that minimizes
## the variance

default.QP <- function(m, dmu=FALSE, mu.constr=NULL)
  {
    ## check the m argument
    if(!is.numeric(m) && length(m) == 1 && m > 0)
      stop("number of cols, m, should be a positive integer scalar")
    
    ## the sum of weights must be equal to 1
    Amat <- matrix(rep(1, m), ncol=1)
    b0 <- 1
    meq <- 1

    ## each w one must be positive
    Amat <- cbind(Amat, diag(rep(1, m)))
    b0 <- c(b0, rep(0, m))

    ## each one less than 1
    Amat <- cbind(Amat, diag(rep(-1, m)))
    b0 <- c(b0, rep(-1, m))

    ## construct the dvec and dmu
    if(!is.logical(dmu) || length(dmu) != 1)
      stop("dmu should be a scalar logical")
    if(dmu) dvec <- rep(1, m)
    else dvec <- rep(0, m)

    ## check for a constraint on mu
    if(!is.null(mu.constr)) {
      if(!is.numeric(mu.constr))
        stop("mu.constr should numeric or NULL")

      ## create columns to add on to Amat that are ones that
      ## alternate in sign 
      addc <- matrix(1, ncol=length(mu.constr), nrow=m)
      parity <- rep(c(1,-1), ceiling(ncol(addc)/2))
      if(ceiling(ncol(addc)/2) != floor(ncol(addc)/2))
        parity <- parity[-length(parity)]
      parity <- matrix(parity, nrow=nrow(Amat), ncol=length(parity), byrow=TRUE)
      addc <- addc * parity

      ## add those columns to Amat, and the constrants to b0
      Amat <- cbind(Amat, addc)
      b0 <- c(b0, mu.constr)

      ## record the length and which colums were added
      mu.constr <- c(length(mu.constr),
                     (ncol(Amat)-(length(mu.constr)-1)):ncol(Amat))
    } else mu.constr <- 0

    ## return the list
    return(list(m=m, dvec=dvec, dmu=dmu, Amat=Amat, b0=b0,
                mu.constr=mu.constr, meq=meq))
  }


## postprocess.QP:
##
## re-arrange the components of QP and the output W
## according to the factors and monomvn ordering used
## for passing back to the calling environment

postprocess.QP <- function(QP, W, M, T, nam)
{
  ## save the QP inputs, and NULL some out
  QP$q <- NULL
  cols <- QP$cols; QP$cols <- NULL
  o <- QP$o
  
  ## convert the W-vector into a matrix
  W[abs(W) < .Machine$double.eps] <- 0
  W <- matrix(W, ncol=length(cols), nrow=T, byrow=TRUE)
  if(!is.null(o)) { ## reorder rows or columns
    oo <- order(o)
    QP$Amat <- QP$Amat[oo,]
    W <- W[,oo,drop=FALSE]
  }

  ## deal with names of columns
  if(! is.null(nam)) {
    if(M-length(cols) > 0) ## remove the factor names
      namr <- nam[-(1:(M-length(cols)))]
    else namr <- nam
    rownames(QP$Amat) <- namr
    colnames(W) <- namr
  }
  
  ## return QP and W as a list
  return(list(QP=QP, W=W))
}


## monomvn.solve.QP:
##
## function to use the QP structure to pass arguments
## top the real solve.QP function

monomvn.solve.QP <- function(S, QP, mu = NULL)
{
  ## check the dims of S
  if(!is.matrix(S) || ncol(S) != nrow(S))
    stop("must have a symmetric matrix for S")

  ## optionally check the dims of mu
  if(!is.null(mu) && length(mu) != ncol(S))
    stop("must have length(mu) = ncol(S)")

  ## check that Amat agrees with dim of ncol(S)
  if(ncol(S) != nrow(QP$Amat))
    stop("must have ncol(S) = nrow(QP$Amat)")
  
  ## check the dimension of S versus QP
  QP <- check.QP(QP, ncol(S), NULL)

  ## apply the augmentations specified by QP

  ## make mu/dvec adjustment
  if(QP$dmu) {
    if(is.null(mu)) stop("QP requires that mu be specified")
    QP$dvec <- QP$dvec * mu
  }

  ## fill rows of Amat with mu
  muc <- QP$mu.constr[1]
  if(muc > 0) {
    mu.constr <- QP$mu.constr[-1]

    ## sanity check 
    if(length(mu.constr) != muc)
      stop("muc v. length(QP$mu.constr[-1]) mismatch")
       
    ## loop over the mu-constraints
    for(i in 1:muc) {
      QP$Amat[,mu.constr[i]] <- QP$Amat[,mu.constr[i]] * mu
    }
  }

  ## finally, solve the Quadratic Program
  sol <- solve.QP(S, QP$dvec, QP$Amat, QP$b0, meq=QP$meq)$sol

  ## set small values of the solultion to zero
  sol[abs(sol) < .Machine$double.eps] <- 0

  ## return the solution
  return(sol)
}


## add.pe.QP:
##
## adds points to an existing plot showing the Bayesian mean
## and MLE quadratic programming solutions for comparison

add.pe.QP <- function(b, m)
  {
    ## sainty check b
    if(is.null(b$QP)) stop("b does not have a QP field")

    ## calculate which initial columns should be removed as factors 
    rem <- ncol(b$S) - b$QP$m
    if(rem < 0) stop("non-conforming objects (b,m) provided")
    cols <- 1:ncol(b$S)
    if(rem > 0) cols <- cols[-(1:rem)]

    ## calculate the mean and mle QP solutions
    W.mean <- monomvn.solve.QP(b$S[cols,cols], b$QP)
    W.mle <- monomvn.solve.QP(m$S[cols,cols], b$QP)

    ## add these points to an existing plot
    points(apply(b$W, 2, mean), col=3, pch=18)
    points(W.mean, col=4, pch=17)
    points(W.mle, col=5, pch=16)

    ## add a legend to the plot
    legend("topright",
           c("MAP", "mean(W)", "w(S.mean)", "MLE"),
           col=2:5, pch=c(21,18:16))

  }
