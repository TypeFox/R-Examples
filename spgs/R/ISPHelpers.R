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

ISPX2vec <- function(x)
#make sure x can be interpreted as a valid stochastic vector of length 4 of the form (A,C,G,T)
{
#Check type of x
	if (!is.vector(x, "numeric"))
		stop("x must be a numeric vector")

#Check element names (if any) of x
	enames <- tolower(names(x)) #names of x
	if (length(enames)>0)
	{
		vec <- c(a=0, c=0, g=0, t=0) #empty, template vector of named relative frequencies
		nucs <- names(vec) #set of nucleic acids
		names(x) <- enames #convert names to lowercase
		common.names <- intersect(nucs, enames) #nucleic acids x and nucs have in common
		if (length(common.names)==0)
			stop("x has named entries, but none corresponding to a, c, g or t")
		vec[common.names] <- x[common.names] #copy relative frequencies of common nucleic acids from x to template
		vec <- vec/sum(vec) #normalise to sum to unity, in case there were extra entries in x
	}
	else
	{
		if (length(x)<4) stop("x contains fewer than 4 entries")
		vec <- x[1:4] #copy x to vec, truncating to 4 elements if necessary
		vec <- vec/sum(vec) #normalise to sum to unity
	} #if
#Check that vec is a non-zero stochastic vector
	if (any(vec<0) || all(vec<=0)) stop("x must be non-negative with at least one positive entry")
	vec
} #function



ISPX2mat <- function(x)
#make sure x can be interpreted as a valid 4X4 stochastic matrix with rows and columns in the order A, C, G, T
{
#Check type and size of x
	if (!is.numeric(x) || !is.matrix(x) || any(dim(x)<4))
		stop("x must be a 4 X 4 numeric matrix/array")

#Check row and column names of x
	nucs <- c("a", "c", "g", "t") #set of nucleic acids
	rnames <- tolower(rownames(x)) #row names of x
	cnames <- tolower(colnames(x)) #column names of x
	if (length(rnames)>0 && length(cnames)>0)
	{
		mat <- matrix(0, 4, 4, dimnames=list(nucs, nucs)) #empty, template vector of named relative dinucleotide frequencies
		common.rnames <- intersect(nucs, rnames) #row names x and mat have in common
		common.cnames <- intersect(nucs, cnames) #column names x and mat have in common
		if (length(common.rnames)==0) stop("x has row names, but none corresponding to a, c, g or t")
		if (length(common.cnames)==0) stop("x has column names, but none corresponding to a, c, g or t")
		#Build a matrix of dinucleotide pairs (1 pair per row) for x
		crEntries <- t(sapply(sapply(as.vector(outer(cnames, rnames, paste0)), function(x) strsplit(x, "")), c))
		mat[crEntries] <- x[crEntries] #copy relative dinucleotide frequencies from x into the template
		mat <- mat/rowSums(mat) #normalise all row sums to unity, in case there were extra entries in x
	} #if
	else
	{
		if (dim(x)[1]<4) stop("x has less than 4 rows")
		if (dim(x)[2]<4) stop("x has less than 4 columns")
		mat <- x[1:4, 1:4] #copy x to mat, truncating to 4X4 if necessary
		mat <- mat/rowSums(mat) #normalise all row sums to unity, in case there were extra entries in x
	} #if
#Check that mat is a non-zero stochastic matrix	
	if (any(mat<0) || all(mat<=0))
		stop("x must be a non-negative matrix with at least one positive entry")
	mat
} #function
