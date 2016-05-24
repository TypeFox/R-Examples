
                                        #******************************************************************************* 
#
# Particle Learning of Gaussian Processes
# Copyright (C) 2010, University of Cambridge
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


## covar:
##
## calculate the covarance matrix betteen the rows of X1
## and X2 with lengthscale d and nugget g

covar <- function(X1=NULL, X2=NULL, d, g)
  {
    ## coerse arguments and extract dimensions
    if(!is.matrix(X1)) {
      if(is.null(X2)) stop("X2 cannot be NULL in this context")
      m <- ncol(X2)
      X1 <- matrix(X1, ncol=m)
    } else m <- ncol(X1)
    n1 <- nrow(X1)

    ## check d and g arguments
    if(length(d) != 1) stop("bad d argument")
    if(length(g) != 1) stop("bad g argument")

    if(is.null(X2)) {

      ## calculate D 
      D <- distance(X1)
      
      ## calculate K
      outK <- .C("dist2covar_symm_R",
                 D = as.double(D),
                 n = as.integer(n1),
                 d = as.double(d),
                 g = as.double(g),
                 K = double(n1 * n1),
                 PACKAGE="plgp")
      
      ## return the covariance matrix
      return(matrix(outK$K, ncol=n1, byrow=TRUE))
      
    } else {
    
      ## check inputs
      if(ncol(X1) != ncol(X2)) stop("col dim mismatch for X1 & X2")
      
      ## coerse arguments and extract dimensions
      X2 <- as.matrix(X2)
      n2 <- nrow(X2)
      
      ## calculate D 
      D <- distance(X1, X2)
      
      ## calculate K
      outK <- .C("dist2covar_R",
                 D = as.double(t(D)),
                 n1 = as.integer(n1),
                 n2 = as.integer(n2),
                 d = as.double(d),
                 g = as.double(g),
                 K = double(n1 * n2),
                 PACKAGE="plgp")
      
      ## return the covariance matrix
      return(matrix(outK$K, ncol=n2, byrow=TRUE))
    }
  }


## covar.sep:
##
## calculate the separable covarance matrix betteen the rows
## of X1 and X2 with lengthscale d and nugget g

covar.sep <- function(X1=NULL, X2=NULL, d, g)
  { 
    ## coerse arguments and extract dimensions
    if(!is.matrix(X1)) {
      if(is.null(X2)) stop("X2 cannot be NULL in this context")
      m <- ncol(X2)
      X1 <- matrix(X1, ncol=m)
    } else m <- ncol(X1)
    n1 <- nrow(X1)

    ## check d and g arguments
    if(length(d) != m) stop("bad d argument")
    if(length(g) != 1) stop("bad g argument")

    if(is.null(X2)) {

      ## calculate K
      outK <- .C("covar_sep_symm_R",
                 col = as.integer(m),
                 X = as.double(t(X1)),
                 n = as.integer(n1),
                 d = as.double(d),
                 g = as.double(g),
                 K = double(n1 * n1),
                 PACKAGE="plgp")
      
      ## return the covariance matrix
      return(matrix(outK$K, ncol=n1, byrow=TRUE))
      
    } else {
    
      ## check inputs
      if(ncol(X1) != ncol(X2)) stop("col dim mismatch for X1 & X2")
      
      ## coerse arguments and extract dimensions
      X2 <- as.matrix(X2)
      n2 <- nrow(X2)
      
      ## calculate K
      outK <- .C("covar_sep_R",
                 col = as.integer(m),
                 X1 = as.double(t(X1)),
                 n1 = as.integer(n1),
                 X2 = as.double(t(X2)),
                 n2 = as.integer(n2),
                 d = as.double(d),
                 g = as.double(g),
                 K = double(n1 * n2),
                 PACKAGE="plgp")
      
      ## return the covariance matrix
      return(matrix(outK$K, ncol=n2, byrow=TRUE))
    }
  }


## covar.sim:
##
## calculate the sim covarance matrix betteen the rows
## of X1 and X2 with index parameter d and nugget g

covar.sim <- function(X1=NULL, X2=NULL, d, g)
  { 
    ## coerse arguments and extract dimensions
    if(!is.matrix(X1)) {
      if(is.null(X2)) stop("X2 cannot be NULL in this context")
      m <- ncol(X2)
      X1 <- matrix(X1, ncol=m)
    } else m <- ncol(X1)
    n1 <- nrow(X1)

    ## check d and g arguments
    if(length(d) != m) stop("bad d argument")
    if(length(g) != 1) stop("bad g argument")

    if(is.null(X2)) {

      ## calculate K
      outK <- .C("covar_sim_symm_R",
                 col = as.integer(m),
                 X = as.double(t(X1)),
                 n = as.integer(n1),
                 d = as.double(d),
                 g = as.double(g),
                 K = double(n1 * n1),
                 PACKAGE="plgp")
      
      ## return the covariance matrix
      return(matrix(outK$K, ncol=n1, byrow=TRUE))
      
    } else {
    
      ## check inputs
      if(ncol(X1) != ncol(X2)) stop("col dim mismatch for X1 & X2")
      
      ## coerse arguments and extract dimensions
      X2 <- as.matrix(X2)
      n2 <- nrow(X2)
      
      ## calculate K
      outK <- .C("covar_sim_R",
                 col = as.integer(m),
                 X1 = as.double(t(X1)),
                 n1 = as.integer(n1),
                 X2 = as.double(t(X2)),
                 n2 = as.integer(n2),
                 d = as.double(d),
                 g = as.double(g),
                 K = double(n1 * n2),
                 PACKAGE="plgp")
      
      ## return the covariance matrix
      return(matrix(outK$K, ncol=n2, byrow=TRUE))
    }
  }


## dist2covar.symm:
##
## given a distance matrix, convert it into a
## a covariance matrix

dist2covar.symm <- function(D, d, g)
  {
    ## extract dimensions and sanity check
    n <- nrow(D)
    if(n != ncol(D)) stop("D must be symmetric")
    
    ## calculate K
    outK <- .C("dist2covar_symm_R",
               D = as.double(t(D)),
               n = as.integer(n),
               d = as.double(d),
               g = as.double(g),
               K = double(n * n),
               PACKAGE="plgp")
    
     ## return the covariance matrix
    return(matrix(outK$K, ncol=n, byrow=TRUE))
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
                 PACKAGE="plgp")
      
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
                 PACKAGE="plgp")
      
      ## return the distance matrix
      return(matrix(outD$D, ncol=n2, byrow=TRUE))
    }
  }
