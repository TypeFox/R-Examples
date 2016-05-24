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


#Extended Pareto distribution
dEPD <- function(x, eta, delta, rho, tau, log = FALSE) #parametrization rho = -tau eta; non zero for x>1
{
    if(!missing(tau) && !missing(rho))
        stop("both tau and rho are defined.")
    if(missing(rho) && !missing(tau))
        rho <- -tau*eta
    if(!missing(rho) && missing(tau))
        tau <- -rho/eta
    if(missing(tau) && missing(rho))
        stop("tau or rho must be defined.")
        
	if(eta <= 0 || delta <= max(-1, eta/rho) || tau <= 0)
        return(rep(NA, length(x)))
	
	dx <- 1/eta * x^(-1/eta-1) * (1 + delta*(1-x^(rho/eta)) )^(-1/eta-1)
	dx <- dx * (1 + delta*(1-(1+rho/eta)*x^(rho/eta)) )
	dx[is.nan(dx)] <- 0
    if(log)
        return( log(dx) )
    else
        return( dx )
}



pEPD <- function(q, eta, delta, rho, tau, lower.tail=TRUE, log.p = FALSE) #parametrization rho = -tau eta; non zero for q>1
{
    if(!missing(tau) && !missing(rho))
        stop("both tau and rho are defined.")
    if(missing(rho) && !missing(tau))
        rho <- -tau*eta
    if(!missing(rho) && missing(tau))
        tau <- -rho/eta
    if(missing(tau) && missing(rho))
        stop("tau or rho must be defined.")
    
	if(eta <= 0 || delta <= max(-1, eta/rho) || tau <= 0)
        return(rep(NA, length(q)))
	
	dq <- 1 - (q * (1 + delta - delta*q^(rho/eta)))^(-1/eta)
	dq[is.nan(dq) | q <= 1] <- 0
	if(!lower.tail)
        dq <- 1-dq
    
    if(log.p)
        return( log(dq) )
    else
        return( dq )
}

#(1-p)^(-eta) = x * (1 + delta - delta*x^(rho/eta)) minimize as a function of x
qEPD <- function(p, eta, delta, rho, tau, lower.tail=TRUE, log.p = FALSE,
    control=list())
{
    if(!missing(tau) && !missing(rho))
        stop("both tau and rho are defined.")
    if(missing(rho) && !missing(tau))
        rho <- -tau*eta
    if(!missing(rho) && missing(tau))
        tau <- -rho/eta
    if(missing(tau) && missing(rho))
        stop("tau or rho must be defined.")
        
    #default control parameters
    con <- list(upperbound=1e6, tol=1e-9)
    namc <- names(con)
    con[namc <- names(control)] <- control

    if(log.p)
        p <- exp(p)
    
    n <- length(p)
    sp <- sort(p, decreasing=TRUE)
    idp <- cbind(1:n, order(p, decreasing=TRUE))
    revidp <- idp[order(idp[, 2]), 1]

    x <- numeric(n)
    x[1] <- getqEPD(sp[1], eta, delta, rho, upperbound=con$upperbound)
    
    if(n >= 2)
    {
        for(i in 2:n)
        x[i] <- getqEPD(sp[i], eta, delta, rho,
            upperbound=ifelse(is.na(x[i-1]), con$upperbound, x[i-1]), tol=con$tol)
    
        x <- x[revidp]
    }
    #null probabilities : lower bound of the domain
    x[abs(p) < .Machine$double.eps] <- 1
    #one probabilities : upper bound of the domain
    x[abs(1-p) < .Machine$double.eps] <- Inf
    x[p < 0 | p > 1] <- NaN
    x
}

getqEPD <- function(p, eta, delta, rho, upperbound=Inf, tol=1e-9)
{
    L2 <- function(x)
    {
        ((1-p)^(-eta) - x * (1 + delta - delta * x^(rho/eta)))^2
    }
    res <- try(optimize(L2, c(1, upperbound), tol=tol), silent=TRUE)
    if(class(res) == "try-error")
        return(NA)
    else
        return(res$minimum)
}


rEPD <- function(n, eta, delta, rho, tau)
{
    NULL
}

#unit pareto distribution
dupareto <- function(x, log=FALSE)
{
    dx <- 1/x^2
    dx[x <= 0] <- 0
    if(log)
        return( log(dx) )
    else
        return( dx )
}

pupareto <- function(q, lower.tail=TRUE, log.p = FALSE)
{
	
	dq <- 1 - 1/q
	dq[q <= 0] <- 0
	if(!lower.tail)
        dq <- 1-dq
    
    if(log.p)
        return( log(dq) )
    else
        return( dq )
}


qupareto <- function(p, lower.tail=TRUE, log.p = FALSE)
{
    if(log.p)
        p <- exp(p)
    if(!lower.tail)
        p <- 1-p
    dq <- 1/(1-p)
    
    #null probabilities : lower bound of the domain
    dq[abs(p) < .Machine$double.eps] <- 0
    #one probabilities : upper bound of the domain
    dq[abs(1-p) < .Machine$double.eps] <- Inf
    dq[p < 0 | p > 1] <- NaN
    
    return(dq)
}

rupareto <- function(n)
    1/runif(n, 0, 1)



#Frechet distribution
dfrechet <- function(x, shape, xmin, log = FALSE)
{
    if(missing(shape))
        stop("missing shape argument")
    if(missing(xmin))
        stop("missing xmin argument")
    if(shape <= 0)
        stop("shape argument must be strictly positive.")
    
	dx <- exp(-(x - xmin)^(-shape))
	dx <- dx * shape * (x - xmin)^(-shape-1)
	dx[is.nan(dx)] <- 0
    if(log)
        return( log(dx) )
    else
        return( dx )
    
}


pfrechet <- function(q, shape, xmin, lower.tail=TRUE, log.p = FALSE)
{
    if(missing(shape))
        stop("missing shape argument")
    if(missing(xmin))
        stop("missing xmin argument")
    if(shape <= 0)
        stop("shape argument must be strictly positive.")
 
	p <- exp(-(q - xmin)^(-shape))
    p[is.nan(p) | q <= xmin] <- 0
    if(!lower.tail)
        p <- 1-p

    if(log.p)
        return( log(p) )
    else
        return( p )
}

qfrechet <- function(p, shape, xmin, lower.tail=TRUE, log.p = FALSE)
{
    if(missing(shape))
        stop("missing shape argument")
    if(missing(xmin))
        stop("missing xmin argument")
    if(shape <= 0)
        stop("shape argument must be strictly positive.")
    if(log.p)
        p <- exp(p)
    if(!lower.tail)
        p <- 1-p
    dq <- xmin+(-log(p))^(-1/shape)

    #null probabilities : lower bound of the domain
    dq[abs(p) < .Machine$double.eps] <- xmin
    #one probabilities : upper bound of the domain
    dq[abs(1-p) < .Machine$double.eps] <- Inf
    dq[p < 0 | p > 1] <- NaN

    return(dq)
}


rfrechet <- function(n, shape, xmin)
{
    if(missing(shape))
        stop("missing shape argument")
    if(missing(xmin))
        stop("missing xmin argument")
    if(shape <= 0)
        stop("shape argument must be strictly positive.")
    if(length(n) > 1)
        n <- length(n)
    u <- runif(n, 0, 1)
    qfrechet(u, shape=shape, xmin=xmin)
}




#unit Frechet distribution
dufrechet <- function(x, log = FALSE)
{
    
	dx <- exp(-1/x) * (x)^(-2)
	dx[is.nan(dx)] <- 0
    if(log)
        return( log(dx) )
    else
        return( dx )
}


pufrechet <- function(q, lower.tail=TRUE, log.p = FALSE)
{
    
	p <- exp(-1/q)
    p[is.nan(p) | q <= 0] <- 0
    if(!lower.tail)
        p <- 1-p
    
    if(log.p)
        return( log(p) )
    else
        return( p )
}

qufrechet <- function(p, lower.tail=TRUE, log.p = FALSE)
{

    if(log.p)
        p <- exp(p)
    if(!lower.tail)
        p <- 1-p
    dq <- 1/(-log(p))
    
    #null probabilities : lower bound of the domain
    dq[abs(p) < .Machine$double.eps] <- 0
    #one probabilities : upper bound of the domain
    dq[abs(1-p) < .Machine$double.eps] <- Inf
    dq[p < 0 | p > 1] <- NaN
    
    return(dq)
}


rufrechet <- function(n)
{

    if(length(n) > 1)
        n <- length(n)
    u <- runif(n, 0, 1)
    qufrechet(u)
}


