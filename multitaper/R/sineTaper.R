##     The multitaper R package
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2011 Karim Rahim 
##
##     Written by Wesley Burr.
##
##     This file is part of the multitaper package for R.
##     http://cran.r-project.org/web/packages/multitaper/index.html
## 
##     The multitaper package is free software: you can redistribute it and 
##     or modify it under the terms of the GNU General Public License as 
##     published by the Free Software Foundation, either version 2 of the 
##     License, or any later version.
##
##     The multitaper package is distributed in the hope that it will be 
##     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
##     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##
##     You should have received a copy of the GNU General Public License
##     along with multitaper.  If not, see <http://www.gnu.org/licenses/>.
##
##     If you wish to report bugs please contact the author:
## 
##     Wesley Burr
##     wesley.burr@gmail.com



####################################################################
##
##  sineTaper
##
##  Generates k sine tapers of length n. These are not actually
##  used in the \pkg{multitaper} implementation of the sine-tapered
##  multiple taper spectrum estimate, but are provided for 
##  plotting and transfer function purposes.
##
##  ref: Kurt S. Riedel and Alexander Sidorenko 
## 
####################################################################

sineTaper <- function(n, k) {

    stopifnot(n >= 8, k >= 1)

    coef1 <- as.double(sqrt(2/(n+1)))
    coef2 <- as.double((pi/(n+1))*seq(1,n,1))
    kmat <- matrix(data=as.double(rep(seq(1,n),each=k)),nrow=n,ncol=k)
    
    taper <- coef1*sin(coef2*kmat)

    out <- NULL
    out$v <- as.matrix(taper) 

    # include k in object since these tapers are not always computed
    # in context of a mtm object.
    res <- list(v=out$v,
                eigen=NULL,
                k=k)
    class(res) <- "dpss"
    return(res)
}
