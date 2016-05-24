#############################################################################
#   Copyright (c) 2014 Christophe Dutang
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the
#   Free Software Foundation, Inc.,
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
#
#############################################################################


#minimum distance power divergence
MDPD <- function(theta, densfun, obs, alpha, ..., control=list())
{
    if(!is.function(densfun))
        stop("densfun should be a density function")
    if(!is.vector(obs) || !is.numeric(obs))
        stop("observations should a numeric vector.")
    if(length(alpha) != 1)
        stop("alpha should be a positive scalar.")
    if(alpha < 0)
        stop("alpha should be a positive scalar.")
    
    
    
    #wrapping function
    f <- function(x)
        do.call(densfun, c(list(x), as.list(theta), ...))
        
    #sanity check
    testf <- try(f(obs))
    if(class(testf) == "try-error")
        stop("the density function raised an error.")
    
    
	n <- length(obs)

    #default control parameters
    con <- list(eps=1e-3, tol=1e-3, lower=1, upper=Inf)
    namc <- names(con)
    con[namc <- names(control)] <- control

	if(alpha == 0)
	{
		res <- -sum( log( f(obs) ) ) / n
	}else
	{
        res <-  -(1 + 1 / alpha) * sum( f(obs)^(alpha) ) / n
		
        g <- function(x)
            do.call(densfun, c(list(x), as.list(theta), ...) )^(alpha+1)
		I <- integrate(g, lower=con$lower, upper=con$upper, rel.tol=con$tol, stop.on.error=FALSE)
        
		if(I$message != "OK") #roundoff error
		{
			print(I)
			print(c(theta, alpha, con$tol, con$eps))
			stop("divergence of integral computation.")
		}else if(I$value < con$tol) #numerical instability
		{
			warning(paste("numerical instability for computing the integral: tol=", I$value))
			I$value <- con$eps * g(con$lower+con$eps) #compute an approximate value
		}
		res <- res + I$value
	}
	res
}

