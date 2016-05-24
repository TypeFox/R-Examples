#This file is part of the source code for
#SPGS: an R package for identifying statistical patterns in genomic sequences.
#Copyright (C) 2015  Universidad de Chile and INRIA-Chile
#
#This program is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 2 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#A copy of Version 2 of the GNU Public License is available in the 
#share/licenses/gpl-2 file in the R installation directory or from 
#http://www.R-project.org/Licenses/GPL-2.

chargaff.gibbs.test <- function(x, maxLag=200)
{
	#Check arguments
	dname <- deparse(substitute(x))
	if (!is.numeric(maxLag) || length(maxLag)!=1 || maxLag!=floor(maxLag) || maxLag<10)
  	stop("maxLag must be an integer greater than or equal to 10")

	if (is.array(x))
	{
		if (length(dim(x))!=5 || !identical(dim(x)[-5], rep(4L, 4)))
			stop("x is an invalid array")
  	nLags <- dim(x)[5]-1
  	maxLag <- nLags-1
  	if (maxLag<10)stop("The fifth dimension of the cylinder count array must be greater than or equal to 10")
		quads <- x
		n <- sum(quads)/nLags
		for (k in 1:nLags)
		{
			quadk <- quads[,,,,k]
			if (sum(quadk)!=n)
				stop("invalid cylinder count array at lag", format(k-1))
		} #for k
#Test cylinder count invariants
		prepairs <- apply(quads[,,,,1], 1:2, sum)
		postpairs <- apply(quads[,,,,1], 3:4, sum)
		if (!identical(prepairs, postpairs)) stop("there is a discrepancy in the counts of 2-cylinders at lag 0")
		for (k in 2:nLags)
		{
			prepairsk <- apply(quads[,,,,k], 1:2, sum)
			postpairsk <- apply(quads[,,,,k], 3:4, sum)
			if (!identical(prepairsk, prepairs) || !identical(postpairsk, postpairs)) stop("There is a discrepancy in the counts of 2-cylinders at lag", format(k-1))
		} #for k
		covInfo <- FHatCovariances(quads)
	}
	else if (is.character(x)) #is x a character vector?
	{ #yes, input is a sequence
#Count cylinders
  lag <- 10
  cutoff <- 0
  while (cutoff==0 || (cutoff>=lag && lag<maxLag))
  {
		lag <- min(c(2*lag, maxLag))
		quads <- cyl2lag2.counts(x, lag)
		covInfo <- FHatCovariances(quads)
    cutoff = covInfo$cutoff
	} #while
	pairs <- apply(quads[,,,,1], 1:2, sum)
	} #if
	else
		stop("x must be a character vector representing a DNA sequence")

#Compute n\hat F vector
	if (!identical(rowSums(pairs), colSums(pairs)))
  	stop("there is a discrepancy in the base counts")
	n <- length(x)
	fhat <- rep(0,5)
	fhat[1] <- pairs[1,1] - pairs[4,4]
	fhat[2] <- pairs[1,2] - pairs[3,4]
	fhat[3] <- pairs[1,3] - pairs[2,4]
	fhat[4] <- pairs[2,1] - pairs[4,3]
	fhat[5] <- pairs[2,2] - pairs[3,3]

#Compute test statistic
	VInv <- try(solve(covInfo$v), silent=TRUE) #get the inverse of the asymptotic variance-covariance matrix
	if (class(VInv)=="try-error")
	{
		stat <- NA
		warning("It is not possible to calculate the test statistic since the covariance matrix v estimated for fhat is singular.")
	}
	else
	{
		stat <- c(t(fhat)%*%VInv%*%fhat/n) #the \chi_5^2 statistic under the null
		if (stat<0)
			warning("The calculated test statistic is negative and meaningless since the covariance matrix v estimated for fhat is not positive definite.")
	} #if
	names(stat) <- "eta"

#Compute p-value
	if (!is.na(stat) && stat>=0)
		p.value <- pchisq(stat, df=5, lower.tail=FALSE)
	else
		p.value <- NA
	names(p.value) <- "p-value"
	
#Return result
	method <- "Test of Chargaff's second parity rule (CSPR) for dinucleotides for Gibbsian sequences\n"
	rval <- list(statistic=stat, p.value=p.value, method=method,
		data.name=dname, FHat=fhat, v=covInfo$v, pairs=pairs, n=n, cutoff=covInfo$cutoff, 
		maxLag=maxLag)
	class(rval) <- "htest"
	rval
} #function
