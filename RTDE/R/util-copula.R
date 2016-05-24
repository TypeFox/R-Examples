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


#FGM copula : Eyraud Farlie Gumbel Morgenstern
dFGM <- function(u, v, alpha, log = FALSE)
{
    if(missing(alpha))
        stop("missing alpha parameter")
    if(alpha < -1 || alpha > 1)
        stop("alpha should be in [-1,1].")

    duv <- 1+alpha+4*alpha*u*v-2*alpha*(u+v)
    duv[is.nan(duv)] <- 0
    if(log)
        duv <- log(duv)
        
    return(duv)
}


pFGM <- function(u, v, alpha, lower.tail=TRUE, log.p = FALSE)
{
    if(missing(alpha))
        stop("missing alpha parameter")
    if(alpha < -1 || alpha > 1)
        stop("alpha should be in [-1,1].")
    
    if(lower.tail)
	{
		p <- u*v*(1+alpha*(1-u)*(1-v))
	}else
	{
		u <- 1-u #P(X > x)
		v <- 1-v #P(Y > y)
		#Equation 2.6.1 of Nelsen (2006)
		p <- u+v-1+(1-u)*(1-v)*(1+alpha*u*v)
	}

    p[is.nan(p)] <- 0
    if(log.p)
        p <- log(p)
    
    return(p)
}

qFGM <- function(p, alpha, lower.tail=TRUE, log.p = FALSE)
{
    NULL
}



rFGM <- function(n, alpha)
{
    if(missing(alpha))
        stop("missing alpha parameter")
    if(alpha < -1 || alpha > 1)
        stop("alpha should be in [-1,1].")

    u <- runif(n)
    w <- runif(n)
    a <- 1 + alpha*(1-2*u)
    b <- sqrt(a^2 - 4*(a-1)*w)
    v <- 2*w/(a+b)
    cbind(u, v)
    
}




#Frank copula
dfrank <- function(u, v, alpha, log = FALSE)
{
    if(missing(alpha))
        stop("missing alpha parameter")
    if(alpha <= 0)
        stop("alpha should be strictly positive.")
    if(length(u) != length(v))
        stop("u, v should have the same length")

    eu <- exp(-alpha*u)
    ev <- exp(-alpha*v)
    duv <- alpha * eu * ev /(eu + ev - exp(-alpha) - eu*ev)
    duv[is.nan(duv)] <- 0
    
    if(log)
        duv <- log(duv)
    
    return(duv)
}



pfrank <- function(u, v, alpha, lower.tail=TRUE, log.p = FALSE)
{
    if(missing(alpha))
        stop("missing alpha parameter")
    if(alpha <= 0)
        stop("alpha should be strictly positive.")
    if(length(u) != length(v))
        stop("u, v should have the same length")


    if(lower.tail)
	{
        p <- (1-exp(-alpha*u)) * (1-exp(-alpha*v)) / (1-exp(-alpha))
        p <- -1/alpha*log(1-p)
	}else
	{
		u <- 1-u #P(X > x)
		v <- 1-v #P(Y > y)
		#Equation 2.6.1 of Nelsen (2006)
        p <- (1-exp(-alpha*u)) * (1-exp(-alpha*v)) / (1-exp(-alpha))
		p <- u+v-1-1/alpha*log(1-p)
	}
    
    p[is.nan(p)] <- 0
    if(log.p)
        p <- log(p)
    
    return(p)
}


qfrank <- function(p, alpha, lower.tail=TRUE, log.p = FALSE)
{
    NULL
}

rfrank <- function(n, alpha)
{
    if(missing(alpha))
        stop("missing alpha parameter")
    if(alpha <= 0)
        stop("alpha should be strictly positive.")
    ## reference: Joe (1997, p.147)
    ## numerical instability taken into account
    u <- runif(n)
    w <- runif(n)
    v <- w*(1-exp(-alpha))/(w+ (1 - w)*exp(-alpha*u))
    v <- -1/alpha*log(1 - v)
    cbind(u,v)
}
