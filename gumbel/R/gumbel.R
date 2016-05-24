#############################################################################
#   Copyright (c) 2008 Anne-Lise Caillat, Christophe Dutang,                                                                #
#   Veronique Larrieu and Triet Nguyen                                                                                                      #
#                                                                                                                                                                         #
#   This program is free software; you can redistribute it and/or modify                                               #
#   it under the terms of the GNU General Public License as published by                                         #
#   the Free Software Foundation; either version 2 of the License, or                                                   #
#   (at your option) any later version.                                                                                                            #
#                                                                                                                                                                         #
#   This program is distributed in the hope that it will be useful,                                                             #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                          #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                 #
#   GNU General Public License for more details.                                                                                    #
#                                                                                                                                                                         #
#   You should have received a copy of the GNU General Public License                                           #
#   along with this program; if not, write to the                                                                                           #
#   Free Software Foundation, Inc.,                                                                                                              #
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                                                             #
#                                                                                                                                                                         #
#############################################################################
### ----- File part of R package gumbel -----
###
### the Gumbel Hougaard Copula is an Archimedean copula
### with generator phi(t)=(-ln(t))^alpha 
###


phigumbel <- function(t, alpha=1) (-log(t))^alpha
invphigumbel <- function(t, alpha=1) exp(-t^(1/alpha))

dgumbel <- function(u, v=NULL, alpha, dim=2, warning = FALSE)
{
	#check args                        
	if(alpha < 1 && !is.numeric(alpha))
		stop("invalid argument : alpha\n")
	
	if(is.null(v))
	{
		powAlpha <- function(x) x^(1/alpha-dim)
		
		#pgumbel will reject wrong dimension
		res <- pgumbel(u, v, alpha, dim)
		
		isArray <- (length(dim(u)) == dim)
		if(isArray) #u is dim-dimensional array
			rangeApply <- 1:(dim-1)
		else #u is a matrix with dim columns
			rangeApply <- 1
		
		res <- res * powAlpha(apply(phigumbel(u,alpha), rangeApply, sum))
		
		res <- res * apply(phigumbel(u,alpha-1)/u, rangeApply, prod)
			
		if(dim == 2)
		{	
			res <- res * ( alpha-1 +(apply(phigumbel(u,alpha), rangeApply, sum))^(1/alpha) )
		}
		if(dim == 3)
		{		
			somme <- (apply(phigumbel(u,alpha), rangeApply, sum))^(1/alpha)
			res <- res * ((2*alpha-1)*(alpha-1) + 3*(alpha-1)*somme + somme^2)			
		}
		if(dim > 3)
			stop("dgumbel not yet implemented for dim > 3")
	}
	else
	{
		x <- cbind(u,v)
		return(dgumbel(x, NULL, alpha, 2)) 
	}
	
	
	if(sum(is.nan(res)) && warning) 
	{
		cat("warning NaN produced in dgumbel :",sum(is.nan(res)),"\n")
		cat("aplha ", alpha, "\n")
	}
	
	return(res)
	
}

pgumbel <- function(u, v=NULL, alpha, dim=2)
{
	#check args                        
	if(alpha < 1 && !is.numeric(alpha))
		stop("invalid argument : alpha\n")
	
	if(!is.null(v))
	{	
		x <- cbind(u,v)
		return(pgumbel(x, NULL, alpha, 2)) 
	}

			
	if(is.array(u) || is.matrix(u))              
	{	                                            
		if(length(dim(u)) == dim)                       
		{
			if(dim(u)[dim] != dim)
				stop("invalid dimension\n")	
		}
		else
		{
			if(dim(u)[length(dim(u))] != dim)
				stop("invalid dimension\n")		
		}
		
		if(length(dim(u)) == dim)			               
			return(invphigumbel(apply(phigumbel(u,alpha), 1:(dim-1), sum), alpha))
		else
			return(invphigumbel(rowSums( phigumbel(u,alpha) ), alpha))
	}
	
	if(is.vector(u))
	{	
		if(length(u) != dim)
			stop("invalid dimension\n")
		return(invphigumbel(sum(phigumbel(u,alpha))))
	}


	stop("invalid argument u\n")
	
}

rgumbel <- function(n, alpha, dim=2, method=1)
{
	#check args
	if(alpha < 1 && !is.numeric(alpha))
		stop("invalid argument : alpha\n")
	
	if(method == 1)
	{
		#generate dim*n exponential random variables
		exprand <- matrix( rexp(dim*n), c(n,dim))
	
		#stable random generation S(1/alpha,0,1,0 ; 0) cf. Nolan(2005) and Chambers(1977)
		unifpirand <- runif(n, 0, pi)
		exprand2 <- rexp(n)
		beta <- 1/alpha
		stablerand <- sin((1 - beta) *unifpirand)^((1 - beta)/beta) * (sin(beta * unifpirand)) / (sin(unifpirand))^(1/beta)
		stablerand <- stablerand /( exprand2^(alpha-1) )
	
		#apply the laplace transform
		unifrand <- invphigumbel(exprand/stablerand, alpha)
		if(sum(is.nan(unifrand))) cat("warning NaN produced in rgumbel\n")
	}
	if(method == 2)
	{
		v2 <- runif(n)
		T <- runif(n)
		#K function
		K <- function(t) t-t*log(t)/alpha
		#numerical inverse of K
		Kinv <- function(x) optimize( function(y) (x-K(y))^2 , c(0,1) )$minimum
		v1 <- sapply(T, Kinv)
		#inverse the joint distribution function
		unifrand <- matrix(0, n,2)
		unifrand[,1] <- invphigumbel( phigumbel(v1, alpha)*v2 , alpha)
		unifrand[,2] <- invphigumbel( phigumbel(v1, alpha)*(1-v2) , alpha)
	}
	
	return(unifrand)
}

gumbel.MBE <- function(x, y, marg = "exp")
{
	n <- length(x)
	alphacop<-1/(1-tau(x,y))
	
	if(marg == "exp")
	{
		lambdachap <- 1/mean(x) * (n-1)/n
		muchap <- 1/mean(y) * (n-1)/n
		return(c(lambdachap,muchap,alphacop))
	}
	if(marg == "gamma")
	{
		lambdaX <- mean(x)/var(x)
		alphaX <- mean(x)^2/var(x)
		lambdaY <- mean(y)/var(y)
		alphaY <- mean(y)^2/var(y)
		return(c(lambdaX, alphaX, lambdaY, alphaY, alphacop))
	}
}

gumbel.EML<-function(x, y, marg = "exp")
{
	if(marg == "exp")
	{
		logL <- function(param)
		{
			lambdaX <- param[1]
			lambdaY <- param[2]
			alphaCop <- param[3]
		
			if(lambdaX > 0 && lambdaY > 0 && alphaCop >= 1)
			{	
				res <- sum( remove.naninf( dexp(x, lambdaX, log=TRUE) ) )
				res <- res + sum( remove.naninf( dexp(y, lambdaY, log=TRUE) ) )
				res <- res + sum( remove.naninf( log( dgumbel( pexp(x, lambdaX), pexp(y, lambdaY), alphaCop) ) ) )
			}
			else 
				res <- -Inf
			return(res)
		}
	
		init <- gumbel.MBE(x,y, "exp")
	
		resMaxLike <- optim(init, logL, control=list(fnscale=-1))
		paramchap <- resMaxLike$par
	}	
	if(marg == "gamma")
	{
		logL <- function(param)
		{
			lambdaX <- param[1]
			alphaX <- param[2]
			lambdaY <- param[3]
			alphaY <- param[4]
			alphaCop <- param[5]
	
			if(lambdaX > 0 && lambdaY > 0 && alphaX > 0 && alphaY > 0 && alphaCop >= 1)
			{	
				res <- sum( remove.naninf( dgamma(x, alphaX, lambdaX, log=TRUE) ) )
				res <- res + sum( remove.naninf( dgamma(y, alphaY, lambdaY, log=TRUE) ) )
				res <- res + sum( remove.naninf( log( dgumbel( pgamma(x, alphaX, lambdaX), pgamma(y, alphaY, lambdaY), alphaCop) ) ) )
			}
			else 
				res <- -Inf
			return(res)
		}
	
		init <- gumbel.MBE(x,y, "gamma")
	
		resMaxLike <- optim(init, logL, control=list(fnscale=-1))
		paramchap <- resMaxLike$par
	}
	return(paramchap)
}


gumbel.IFM <- function(x,y, marg = "exp")
{
	if(marg == "exp")
	{
		partLogL<-function(theta, data)
		{
			rate <- theta
			if(rate > 0)
				res <- sum( remove.naninf( log( dexp(data, rate) ) ) )
			else
				res <- -Inf	
			return(res)	
		}		
	
		init<-gumbel.MBE(x, y, "exp")
		range <- c(max(0, init[1] * .5), init[1] * 1.5)
		lambdachap <-optimize(partLogL, range, data=x, maximum=TRUE)$maximum
	
		range <- c(max(0, init[2] * .5), init[2] * 1.5)
		muchap <- optimize(partLogL, range, data=y, maximum=TRUE)$maximum
	
		pseudou <- pexp(x, lambdachap)
		pseudov <- pexp(y, muchap)
	
		copLogL <- function(alpha)
		{
			if(alpha >= 1)
				res <- sum( remove.naninf( log( dgumbel( pseudou, pseudov, alpha) ) ) )
			else
				res <- -Inf	
			return(res)	
		}
	
		alphachap <- optimize(copLogL, c(1+sqrt(.Machine$double.eps),10), maximum=TRUE)$maximum 	
		return(c(lambdachap, muchap, alphachap))			
	}
	
	if(marg == "gamma")
	{
		partLogL<-function(theta, data)
		{
			rate <- theta[1]
			shape <- theta[2]		
			if(rate > 0 && shape > 0)
				res <- sum( remove.naninf( log( dgamma(data, shape, rate) ) ) )
			else
				res <- -Inf	
			return(res)	
		}		
	
		init<-gumbel.MBE(x,y, "gamma")
	
		resMaxLike <- optim(init[1:2], partLogL, control=list(fnscale=-1), data=x)
		thetachap <- resMaxLike$par
		lambdaXchap <- thetachap[1]
		alphaXchap <- thetachap[2]
	
		resMaxLike <- optim(init[3:4], partLogL, control=list(fnscale=-1), data=y)
		thetachap <- resMaxLike$par
		lambdaYchap <- thetachap[1]
		alphaYchap <- thetachap[2]
	
		pseudou <- pgamma(x, alphaXchap, lambdaXchap)
		pseudov <- pgamma(y, alphaYchap, lambdaYchap)
	
		copLogL <- function(alpha)
		{
			if(alpha >= 1)
				res <- sum( remove.naninf( log( dgumbel( pseudou, pseudov, alpha) ) ) )
			else
				res <- -Inf	
			return(res)	
		}
	
		alphacopchap <- optimize(copLogL, c(1+sqrt(.Machine$double.eps),10), maximum=TRUE)$maximum 	
	
		return(c(lambdaXchap, alphaXchap, lambdaYchap, alphaYchap, alphacopchap))			
	}
}

gumbel.CML <- function(x, y)
{
	#pseudo data
	margX <- ecdf(x)
	margY <- ecdf(y)
	
	pseudou <- margX(x)
	pseudov <- margY(y)
	
	copLogL <- function(alpha)
	{
		if(alpha >= 1)
			res <- sum( remove.naninf( log( dgumbel(pseudou, pseudov, alpha) ) ) )
		else
			res <- -Inf	
		
		return(res)	
	}		
	
	return( optimize(copLogL, c(1+sqrt(.Machine$double.eps), 10), maximum=TRUE)$maximum )	
}

#auxiliary functions
remove.naninf <- function(x)
{
	#remove NaN values
	temp <- x[which(!is.nan(x), arr.ind=TRUE)]
	#remove Inf values
	temp[which(is.finite(temp), arr.ind=TRUE)]
}

tau<-function(x,y) return(cor(x,y,method="kendall"))