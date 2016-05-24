#*******************************************************************************
#
# Local Approximate Gaussian Process Regression
# Copyright (C) 2013, The University of Chicago
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
# Questions? Contact Robert B. Gramacy (rbgramacy@chicagobooth.edu)
#
#*******************************************************************************



## rbetter_R:
##
## sample uniformly
## from the (diagonal) half of a rectable where the sum of
## the rows of X (the sample) are less than ystar

rbetter <- function(n, rect, ybest)
  {
    rr <- .C("rbetter_R",
             n = as.integer(n),
             m = as.integer(nrow(rect)),
             rect = as.double(rect),
             ystar = as.double(ybest),
             X = double(n*nrow(rect)),
             PACKAGE = "laGP")

    return(matrix(rr$X, ncol=nrow(rect), byrow=TRUE))
  }


## distance:
##
## calculate the distance matrix between the rows of X1
## and X2, or between x1 and itself when X2=NULL

distance <- function(X1, X2=NULL)
  {
    ## coerse arguments and extract dimensions
    X1 <- as.matrix(X1)
    n1 <- nrow(X1)
    m <- ncol(X1)

    if(is.null(X2)) {

      ## calculate D
      outD <- .C("distance_symm_R",
                 X = as.double(t(X1)),
                 n = as.integer(n1),
                 m = as.integer(m),
                 D = double(n1 * n1),
                 PACKAGE = "laGP")
      
      ## return the distance matrix
      return(matrix(outD$D, ncol=n1, byrow=TRUE))
      
    } else {

      ## coerse arguments and extract dimensions
      X2 <- as.matrix(X2)
      n2 <- nrow(X2)

      ## check inputs
      if(ncol(X1) != ncol(X2)) stop("col dim mismatch for X1 & X2")
      
      ## calculate D
      outD <- .C("distance_R",
                 X1 = as.double(t(X1)),
                 n1 = as.integer(n1),
                 X2 = as.double(t(X2)),
                 n2 = as.integer(n2),
                 m = as.integer(m),
                 D = double(n1 * n2),
                 PACKAGE = "laGP")
      
      ## return the distance matrix
      return(matrix(outD$D, ncol=n2, byrow=TRUE))
    }
  }


## llikGP.d:
##
## build a GP with X,Z,g,d and return the 
## log likelihood, and the derivatives 

llikGP.d <- function(dg, X=NULL, y=NULL, gpi.in=NULL, 
                     ab=c(0,0), param=c("d", "g"))
  {
    ## return error for negative values
    ## if(any(dg <= 1e-6)) return(-Inf)

    ## possibly init GP
    if(is.null(gpi.in)) gpi <- newGP(X,y,dg[1,1],dg[1,2], dK=TRUE)
    else {
      gpi <- gpi.in
      newparamsGP(gpi, dg[1,1], dg[1,2])
      if(!(is.null(X) && is.null(y))) stop("must have gpi or X & y")
    }

    ## check params and reconcile dab and gab
    param <- match.arg(param)
    dab <- gab <- c(0,0)
    if(param == "d") dab <- ab
    else gab <- ab

    ## make llik and dllik calculations
    llik <- rep(NA, nrow(dg))
    dllik <- data.frame(d=llik, d2=llik)
    llik[1] <- llikGP(gpi, dab, gab)
    dllik[1,] <- dllikGP(gpi, ab, param=param)
    if(nrow(dg) > 1) {
      for(i in 2:nrow(dg)) {
        newparamsGP(gpi, dg[i,1], dg[i,2])
        llik[i] <- llikGP(gpi, dab, gab)
        dllik[i,] <- dllikGP(gpi, ab, param=param)
      }
    }

    ## clean up
    if(is.null(gpi.in)) deleteGP(gpi)

    ## return
    return(data.frame(llik, dllik))
  }
