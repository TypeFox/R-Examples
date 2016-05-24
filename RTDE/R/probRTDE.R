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

prob <- function(object, q, ...)
    UseMethod("prob")

prob.default <- function(object, q, ...)
    return(object)

prob.RTDE <- function(object, q, ...)
{
    stopifnot(class(object) == "RTDE")
    
    object$prob <- pRTDE(q, object)
    object$q <- q
    object
}

prob.fitRTDE <- function(object, q, ...)
{
    stopifnot(class(object) == "fitRTDE")
    
    stop("not implemented.")
    
    #object$prob <- pRTDE(q, object)
    #object$q <- q
    #class(object) <- "RTDE"
    #object
}

pRTDE <- function(q, object)
{
    stopifnot(class(object) %in% c("RTDE", "fitRTDE"))
    
    if(class(object) == "RTDE")
    {
        if(length(object$simu) > 0)
            nbreplic <- object$simu$replicate
        else
            nbreplic <- 1
        
        object <- object$fit
    }else
    {
        nbreplic <- 1
    }

    m <- object$m
    n <- object$n
    
    a.length <- length(object$alpha)
    o.length <- length(object$omega)
    m.length <- length(m)
    
    Ztilde <- object$Ztilde
    eta <- object$eta
    delta <- object$delta
    rho <- object$rho
    dnames <- dimnames(eta)
    initletter <- substr(unlist(dnames), 1, 1)
    
    prob <- array(0, dim=c(a.length, o.length, m.length, nbreplic, length(q)))
    eta <- array(eta, dim=c(a.length, o.length, m.length, nbreplic))
    delta <- array(delta, dim=c(a.length, o.length, m.length, nbreplic))
    Ztilde <- array(Ztilde, dim=c(n, o.length, nbreplic))
    
    for(i in 1:a.length)
    {
        for(j in 1:o.length)
        {
            
            for(k in 1:m.length)
            {
                for(l in 1:nbreplic)
                {
                    Ztildej <- Ztilde[, j, l]

                p_mn <- pEPD(q / Ztildej[n-m[k]], eta[i, j, k, l], delta[i, j, k, l], rho=rho, lower.tail=FALSE)
                prob[i, j, k, l, ] <- p_mn * m[k] / n
                
                }
            }
        }
    }
    dimnames(prob) <- list(
        paste("a=", round(object$alpha, 3), sep=""),
        paste("o=", round(object$omega, 3), sep=""),
        paste("m=", m, sep=""),
        paste("s=", 1:nbreplic, sep=""),
        paste("q=", round(q, 3), sep=""))
    



    if(a.length == 1)
    {
        prob <- prob[1,,,,]
    }else if(o.length == 1)
    {
        prob <- prob[,1,,,]
    }else if(m.length == 1)
    {
        prob <- prob[,,1,,]
    }else if(nbreplic == 1)
    {
        prob <- prob[,,,1,]
    }
   
    initletter2 <- substr(unlist(dimnames(prob)), 1, 1)

    if(any(initletter2[initletter2 %in% c("a", "o", "m", "s")] != initletter))
    {
        print(unlist(dnames))
        print(unlist(dimnames(prob)))
        stop("wrong dimnames")
    }
    
    prob
}
