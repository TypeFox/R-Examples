pettitt.test <- function(x){
##    Copyright (C) 2015, 2016  Thorsten Pohlert
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##    This function computes Pettitt test. It calls the Fortran subroutine
##    F_pettitt
##
    na.fail(x)
    n <- length(x)
   
    DNAME <- deparse(substitute(x))
    retval <- list(nobs = n, 
                   statistic = NULL, estimate=NULL,
                   p.value = NULL,
                   data.name=DNAME,
                   alternative="true change point is present in the series",
                   method = "Pettitt's test for single change-point detection")
    res <- .Fortran(F_pettitt, n=as.integer(n),
                x=as.single(x), pval=as.single(0),
                tau=as.integer(0), Kt=as.integer(0))              
    K <- res$Kt
    tau <- res$tau
    pval <- res$pval
    names(K) <- "K"
    retval$statistic <- K
    names(tau) <- "probable change point at tau"
    retval$estimate <- tau
    retval$p.value <-  pval
    class(retval) <- "htest"
    return(retval)
}
