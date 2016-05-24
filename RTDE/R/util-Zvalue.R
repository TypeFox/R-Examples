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


zvalueRTDE <- function(obs, omega, nbpoint,
    output=c("orig", "relexcess"), marg=c("upareto", "ufrechet", "uunif"))
{
    if(!is.matrix(obs) || NCOL(obs) != 2)
        stop("observations should a two-column matrix.")
    if(length(omega) != 1)
        stop("omega should be a scalar.")
    if(missing(nbpoint))
        nbpoint <- 1
    if(length(nbpoint) != 1)
        stop("nbpoint should be an integer.")

    n <- NROW(obs)
    m <- nbpoint
        
    if(m <= 0 || m >= n)
        stop("nbpoint should be an integer in [1, n-1].")
        
    output <- match.arg(output, c("orig", "relexcess"))
    marg <- match.arg(marg, c("upareto", "ufrechet", "uunif"))
    
    rX <- rank(obs[,1])
    rY <- rank(obs[,2])
    
    if(marg == "uunif")
    {
		X <- rX / (n+1)
		Y <- rY / (n+1)
    }else if(marg == "upareto")
    {
        X <- (n+1)/(n+1 - rX)
		Y <- (n+1)/(n+1 - rY)
    }else if(marg == "ufrechet")
    {
        X <- 1/(log((n+1)/rX))
		Y <- 1/(log((n+1)/rY))
    }else
        stop("wrong marginal")
    

    Ztilde <- pmin(X, omega/(1-omega)*Y)
    if(output == "orig")
    {
        res <- list(type="orig", omega=omega, Ztilde = Ztilde, n=n)
    }
    if(output == "relexcess")
    {
        
        Zn_m <- sort(Ztilde)[n-m]
        Z <- Ztilde[Ztilde > Zn_m] / Zn_m
        res <- list(type="relexcess", omega=omega, Z = Z, n=n, m=m)
    }
    class(res) <- "zvalueRTDE"
    res
}



#generic functions
print.zvalueRTDE <- function(x, ...)
{
    if (!inherits(x, "zvalueRTDE"))
        stop("Use only with 'zvalueRTDE' objects")
    
    if(x$type == "orig")
    {
        cat("Z tilde omega variables \n")
        print(x$Ztilde)
    }
    if(x$type == "relexcess")
    {
        cat("Z variables \n")
        print(x$Z)
    }
    
    invisible(x)
}

summary.zvalueRTDE <- function(object, ...)
{
    if (!inherits(object, "zvalueRTDE"))
        stop("Use only with 'zvalueRTDE' objects")
    
    class(object) <- c("summary.zvalueRTDE", class(object))
    object
}

print.summary.zvalueRTDE <- function(x, ...)
{
    if (!inherits(x, "summary.zvalueRTDE"))
        stop("Use only with 'summary.zvalueRTDE' objects")
    
    if(x$type == "orig")
    {
        cat("Z tilde omega variables \n")
        print(summary(x$Ztilde))
    }
    if(x$type == "relexcess")
    {
        cat("Z variables \n")
        print(summary(x$Z))
    }
    
    invisible(x)
}


relexcess <- function(x, nbpoint, ...)
    UseMethod("relexcess")

relexcess.default <- function(x, nbpoint, ...)
    return(x)

relexcess.zvalueRTDE <- function(x, nbpoint, ...)
{
    if (!inherits(x, "zvalueRTDE"))
        stop("Use only with 'zvalueRTDE' objects")
    
    if(missing(nbpoint))
        stop("missing nbpoint argument.")
    if(length(nbpoint) != 1)
        stop("nbpoint should be an integer.")
    if(nbpoint <= 0 || nbpoint >= x$n)
        stop("nbpoint should be an integer in [1, n-1].")
    
    if(x$type == "orig")
    {
        m <- nbpoint
        Zn_m <- sort(x$Ztilde)[x$n-m]
        Z <- x$Ztilde[x$Ztilde > Zn_m] / Zn_m
        
        res <- list(type="relexcess", omega=x$omega, Z = Z, n=x$n, m=m)
        class(res) <- "zvalueRTDE"
        return(res)
    }else
        return(x)
    
}

