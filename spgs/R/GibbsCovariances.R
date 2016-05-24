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

#GibbsCovariances
#Computes the asymptotic covariances between pairs of nucleic acid bases in 
#a genome sequence or using a precomputed  array of counts.
#
#USAGE
#GibbsCovariances(seq, maxLag, epsilon=0.01)
#covs = gibbscovariances0(cylcounts)
#
#INPUTS
#x:  A DNA sequence expressed as a character vector 
#    or an object with class SeqFastadna from the seqinr package.
#    Alternatively, x may be a 4 X 4 X 4 X 4 X (maxLag+1) array of precomputed nucleic 
#    acid counts from a genome.  cylcounts(a, b, c, d, k) contains the 
#    number of times the pairs (a, b) and (c,d) are seen in a genome 
#    seperated by the distance k.  Note that for k=0, (a,b) and (c,d) 
#    overlap exactly so that a would have to equal c and b would have to 
#   equal d in order for cylcounts(a,b,c,d,k) to hold a non-zero value.
#maxLag:  The lag at which to stop summation in the infinite series 
#    representation of the asymptotic covariance.  This parameter is 
#    optional and ddefaults to 10 if it is not specified.
#epsilon:  The factor used to determine at which lag to cut off the 
#    summations used to calculate the covariance matrix between (a,b) 
#    and (c,d).  The default value is 0.01.
#
#OUTPUT
#covs:  A 4 X 4 X 4 X 4 array of estimated asymptotic covariances.
#    covs(a,b,c,d) gives the covariance between the pairs of nucleic acids
#    (a,b) and (c,d).
#
#ALSO SEE
#GibbsCovariances2
#
#DEPENDENCIES
#cyl2lag2.counts

GibbsCovariances <- function(x, maxLag=10, epsilon=0.01)
{
	if (!is.numeric(maxLag) || length(maxLag)>1 || maxLag!=floor(maxLag) || maxLag<10)
	stop("maxLag must be an integer greater than or equal to 10")

	if (is.array(x))
	{
		if (length(dim(x))!=5 || !identical(dim(x)[-5], rep(4L, 4)))
			stop("x is an invalid array")
  	maxLag<- dim(x)[5]
  	if (maxLag<10)stop("The fifth dimension of the cylinder count array must be greater than or equal to 10")
		quads <- x
		n <- sum(quads)/maxLag
	}
	else if (is.character(x))
	{
#Count cylinders
		n <- length(x)
		quads <- cyl2lag2.counts(x, maxLag)
	}
	else
		stop("x must be either a character vector or a 5-dimensional array containing ab.{k}cd counts")

#Compute probabilities
	quadProbs <- quads/n #convert counts to empirical probabilities

#Compute differences of probabilities
	diffs <- quadProbs
	for (a in 1:4)
		for (b in 1:4)
			for (c in 1:4)
				for (d in 1:4)
					diffs[a,b,c,d,] <- diffs[a,b,c,d,] - quadProbs[a,b,a,b,1]*quadProbs[c,d,c,d,1]

	cutoffs <- matrix(0, nrow=4, ncol=4)
	for (a in 1:4)
		for (b in 1:4)
		{
			s <- abs(diffs[a,b,a,b,])
			small <- which(s<=epsilon*s[1])
			cutoffs[a,b] <- ifelse(length(small)>0, min(small),maxLag)
		} #for b,a
	cutoff <- max(c(3, cutoffs))
	diffs <- apply(diffs[,,,,1:cutoff], 1:4, sum)

#Compute the asymptotic covariances \lim_{n \to\infty} Cov(\frac{N_n^{ab}}{\sqrt{n}}, \frac{N_n^{cd}}{\sqrt{n}})
	covs <- array(0, dim=rep(4L, 4))
  	for (a in 1:4)
    	for (b in 1:4)
      	for (c in 1:4)
        	for (d in 1:4)
          	covs[a,b,c,d] <- diffs[a,b,c,d] + diffs[c,d,a,b] - (quadProbs[a,b,c,d,1] - quadProbs[a,b,a,b,1]*quadProbs[c,d,c,d,1])
	list(covs=covs, cutoff=cutoff)
} #function
