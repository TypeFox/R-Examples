##     The spsmooth R package
##     An extension for mgcv and gam.
##     Copyright (C) 2012 Wesley Burr
##
##     This file is part of the spsmooth package for R.
##
##     Written by Wesley Burr and Karim Rahim, taken from 'multitaper'
##     package written by the same.
##
##     The spsmooth package is free software: you can redistribute it 
##     and/or modify it under the terms of the GNU General Public License as 
##     published by the Free Software Foundation, either version 2 of the 
##     License, or any later version.
##
##     The spsmooth package is distributed in the hope that it will be 
##     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
##     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##
##     You should have received a copy of the GNU General Public License
##     along with spsmooth.  If not, see <http://www.gnu.org/licenses/>.
##
##     If you wish to report bugs please contact the author. 
##     Wesley Burr <wburr@mast.queensu.ca>
##     239 Jeffery Hall, Queen's University, Kingston Ontario
##     Canada, K7L 3N6

##########################################################
##
##  .dpss
## 
##  Generates k orthogonal discrete prolate spheroidal
##  sequences (dpss) using the tridiagonal method. See
##  Slepian (1978), page 1379 and Percival and Walden
##  chapter 8.4
##
##########################################################

.dpss <- function(n, k, nw, returnEigenvalues=TRUE) {

    stopifnot(n >= 1, nw >= 0.5, k >= 1, nw/n <= 0.5)
    if(k >= 1.5+2*nw) { warning("Inner-Band concentration compromised by choice of k.") }

    # If k is passed in as floating point, the cast to 
    # as.integer() in the Fortran call does not quite work properly
    if(!is.integer(k)) {
      k<-as.integer(floor(k));
    } 

    ## 'eigen' is of proper length for use by LAPACK functoins.
    ## This will use lapack functions in place of the
    ## EISPACK functions referenced in Percival and Waldern
    out <- .Fortran("dpss", as.integer(n), as.integer(k),
              as.double(nw), 
              v=double(n*k), eigen=double(k),
              PACKAGE='spsmooth')
    out$v <- matrix(data=out$v, nrow=n, ncol=k, byrow=FALSE)
    if(returnEigenvalues) {
        out$eigen <- .dpssToEigenvalues(out$v, nw)
    } else {
        out$eigen <- NULL
    }
    
    res <- list(v=out$v,
                eigen=out$eigen)
    class(res) <- "dpss"
    return(res)
}

##########################################################
##
##  .dpssToEigenvalues
## 
##  Given a set of dpss tapers, find the eigenvalues corresponding
## to the generated dpss's
##
##  See Percival and Walden exercise 8.4
##
##########################################################


.dpssToEigenvalues <- function(v, nw) {
    v <- as.matrix(v)
    n <- length(v[,1])
    k <- length(v[1,])
    w <- nw/n
    npot <- 2**(ceiling(log2(2*n)))

    ## pad
    scratch <- rbind(v, matrix(data=0, nrow=npot-n, ncol=k))
    ## n * acvs
    scratch <- Re(mvfft(abs(mvfft(scratch))**2))/npot
    
    j <- 1:(n-1)
    vectorOfRatios <- c(2*w, sin(2*pi*w*j)/(pi*j))
    
    ##  Note: both vector of ratios and scratch
    ##  roughly decrease in amplitude with increasing index,
    ##  so sum things up in reverse order
    eigenvalues <- NULL
    if(k>1) {
        eigenvalues <- apply(scratch[n:2,] * vectorOfRatios[n:2], 2, sum)
    } else {
        eigenvalues <- sum(scratch[n:2,] * vectorOfRatios[n:2])
    }
    return(2*eigenvalues +
           vectorOfRatios[1]*scratch[1,1:k])
}

