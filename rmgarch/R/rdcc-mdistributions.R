#################################################################################
##
##   R package rmgarch by Alexios Ghalanos Copyright (C) 2008-2013.
##   This file is part of the R package rmgarch.
##
##   The R package rmgarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rmgarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################


.mDistributionBounds = function(distribution)
{
	if (distribution == "mvnorm"){
		skew 	= 0
		skew.LB = 0
		skew.UB = 0
		shape 	= 0
		shape.LB = 0
		shape.UB = 0}
	if (distribution == "mvt"){
		skew 	= 0
		skew.LB = 0
		skew.UB = 10
		shape 	= 5
		shape.LB = 2.01
		shape.UB = 50}
	if (distribution == "mvlaplace"){
		skew 	= 0
		skew.LB = 0
		skew.UB = 0
		shape 	= 1
		shape.LB = 0
		shape.UB = 100}
	
	skewed.dists = NULL
	shaped.dists = c("mvt", "mvlaplace")
	if(any(skewed.dists == distribution)) include.skew = TRUE else include.skew = FALSE
	if(any(shaped.dists == distribution)) include.shape = TRUE else include.shape = FALSE
	return(list(shape = shape, shape.LB = shape.LB, shape.UB = shape.UB, skew = skew,
					skew.LB = skew.LB, skew.UB = skew.UB, include.skew = include.skew, 
					include.shape = include.shape))
}

wmargin = function(distribution = "mvnorm", weights, mean, Sigma, shape = NA, skew = NA){
	ans = switch(distribution,
				mvnorm = wmargin.mvn(weights, Sigma, mean),
				mvlaplace = wmargin.mvlaplace(weights, Sigma, mean),
				mvt = wmargin.mvt(weights, Sigma, mean, shape))
		return( ans )
}
wmargin.mvn = function(weights, cov, mu)
{
	if( !is.array(cov) ) stop("\ncov must be an array")
	n = dim(cov)[3]
	m = dim(cov)[1]
	if( is.matrix( weights ) ){
		mw = dim(weights)[2]
		if( mw != m ) stop("\nInconsistent weights and cov asset dimensions.")
		nw = dim(weights)[1]
		if( nw == 1 ) weights = matrix(weights, ncol = m, nrow = n, byrow = TRUE)
		if( nw != n ) stop("\nInconsistent weights and cov length.")		
	} else{
		mw = length(as.numeric(weights))
		if( mw != m ) stop("\nInconsistent weights and cov asset dimensions.")
		weights = matrix( as.numeric(weights), nrow = n, ncol = m, byrow = TRUE)
	}
		
	if( is.matrix( mu ) ){
		mw = dim(mu)[2]
		if( mw != m ) stop("\nInconsistent mu and cov asset dimensions.")
		nw = dim(mu)[1]
		if( nw == 1 ) mu = matrix(mu, ncol = m, nrow = n, byrow = TRUE)
		if( nw != n ) stop("\nInconsistent mu and cov length.")		
	} else{
		mw = length(as.numeric(mu))
		if( mw != m ) stop("\nInconsistent mu and cov asset dimensions.")
		mu = matrix( as.numeric(mu), nrow = n, ncol = m, byrow = TRUE)
	}
	psigma = sqrt( apply(as.data.frame(1:n), 1, FUN = function(i) weights[i,] %*% cov[,,i] %*% weights[i,]   ) )	
	pmu	= rowSums(weights * mu)
	port = cbind(pmu, psigma, rep(0, n), rep(0, n))
	colnames(port) = c("mu", "sigma", "skew", "shape")
	return( port )
}

wmargin.mvlaplace = function(weights, cov, mu)
{
	if( !is.array(cov) ) stop("\ncov must be an array")
	n = dim(cov)[3]
	m = dim(cov)[1]
	if( is.matrix( weights ) ){
		mw = dim(weights)[2]
		if( mw != m ) stop("\nInconsistent weights and cov asset dimensions.")
		nw = dim(weights)[1]
		if( nw == 1 ) weights = matrix(weights, ncol = m, nrow = n, byrow = TRUE)
		if( nw != n ) stop("\nInconsistent weights and cov length.")		
	} else{
		mw = length(as.numeric(weights))
		if( mw != m ) stop("\nInconsistent weights and cov asset dimensions.")
		weights = matrix( as.numeric(weights), nrow = n, ncol = m, byrow = TRUE)
	}
	
	if( is.matrix( mu ) ){
		mw = dim(mu)[2]
		if( mw != m ) stop("\nInconsistent mu and cov asset dimensions.")
		nw = dim(mu)[1]
		if( nw == 1 ) mu = matrix(mu, ncol = m, nrow = n, byrow = TRUE)
		if( nw != n ) stop("\nInconsistent mu and cov length.")		
	} else{
		mw = length(as.numeric(mu))
		if( mw != m ) stop("\nInconsistent mu and cov asset dimensions.")
		mu = matrix( as.numeric(mu), nrow = n, ncol = m, byrow = TRUE)
	}
	psigma = sqrt( apply(as.data.frame(1:n), 1, FUN = function(i) weights[i,] %*% cov[,,i] %*% weights[i,]   ) )	
	pmu	= rowSums(weights * mu)
	port = cbind(pmu, psigma, rep(0, n), rep(1, n))
	colnames(port) = c("mu", "sigma", "skew", "shape")
	
	return( port )
}

wmargin.mvt = function(weights, cov, mu, dof)
{
	if( !is.array(cov) ) stop("\ncov must be an array")
	n = dim(cov)[3]
	m = dim(cov)[1]
	if( is.matrix( weights ) ){
		mw = dim(weights)[2]
		if( mw != m ) stop("\nInconsistent weights and cov asset dimensions.")
		nw = dim(weights)[1]
		if( nw == 1 ) weights = matrix(weights, ncol = m, nrow = n, byrow = TRUE)
		if( nw != n ) stop("\nInconsistent weights and cov length.")		
	} else{
		mw = length(as.numeric(weights))
		if( mw != m ) stop("\nInconsistent weights and cov asset dimensions.")
		weights = matrix( as.numeric(weights), nrow = n, ncol = m, byrow = TRUE)
	}
	
	if( is.matrix( mu ) ){
		mw = dim(mu)[2]
		if( mw != m ) stop("\nInconsistent mu and cov asset dimensions.")
		nw = dim(mu)[1]
		if( nw == 1 ) mu = matrix(mu, ncol = m, nrow = n, byrow = TRUE)
		if( nw != n ) stop("\nInconsistent mu and cov length.")		
	} else{
		mw = length(as.numeric(mu))
		if( mw != m ) stop("\nInconsistent mu and cov asset dimensions.")
		mu = matrix( as.numeric(mu), nrow = n, ncol = m, byrow = TRUE)
	}
	psigma = sqrt( apply(as.data.frame(1:n), 1, FUN = function(i) weights[i,] %*% cov[,,i] %*% weights[i,]   ) )	
	pmu	= rowSums(weights * mu)
	port = cbind(pmu, psigma, rep(0, n), rep(dof, n) )
	colnames(port) = c("mu", "sigma", "skew", "shape")
	return( port )
}
# use: mnormt
# 
dmdist = function(distribution = "mvnorm", y, mu, sigma, lambda = -0.5, skew, shape) 
{
	valid.distributions = c("mvnorm", "mvt", "mvlaplace", "mvnig", "mvghyp")
	if (!any(valid.distributions == distribution)) stop("\nnot a valid distributions\n", call. = FALSE)
	
}

pmdist = function(distribution = "mvnorm", y, mu, sigma, lambda = -0.5, skew, shape) 
{
	valid.distributions = c("mvnorm", "mvt", "mvlaplace", "mvnig", "mvghyp")
	if (!any(valid.distributions == distribution)) stop("\nnot a valid distributions\n", call. = FALSE)
	
}

qmdist = function(distribution = "mvnorm", y, mu, sigma, lambda = -0.5, skew, shape) 
{
	valid.distributions = c("mvnorm", "mvt", "mvlaplace", "mvnig", "mvghyp")
	if (!any(valid.distributions == distribution)) stop("\nnot a valid distributions\n", call. = FALSE)
	
}

rmdist = function(distribution = "mvnorm", y, mu, sigma, lambda = -0.5, skew, shape) 
{
	valid.distributions = c("mvnorm", "mvt", "mvlaplace", "mvnig", "mvghyp")
	if (!any(valid.distributions == distribution)) stop("\nnot a valid distributions\n", call. = FALSE)
	
}