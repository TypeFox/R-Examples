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

#FHatCovariances
#[Description to be included.]
#
#USAGE
#v = FHatCovariances(cylcounts)
#
#INPUTS
#cylcounts:  A 4X4X4X4XK array of ab.[k-2]cd cylinder counts.
#
#OUTPUT
#The 5X5 variance-covariance matrix for the \hat F vector.
#
#DEPENDENCIES
#GibbsCovariances
#
#ALSO SEE
#etastat
FHatCovariances <- function(cylcounts)
{
	if (length(dim(cylcounts))==4)
	{
		s <- cylcounts
		cutoff <- NA
	}
	else
	{
 		covInfo <- GibbsCovariances(cylcounts)
 		s <- covInfo$covs
 		cutoff <- covInfo$cutoff
	} #if
	v <- matrix(0, 5, 5)
	v[1,2] <- s[1,1,1,2] + s[4,4,3,4] - s[1,1,3,4] - s[4,4,1,2]
	v[1,3] <- s[1,1,1,3] + s[4,4,2,4] - s[1,1,2,4] - s[4,4,1,3]
	v[1,4] <- s[1,1,2,1] + s[4,4,4,3] - s[1,1,4,3] - s[4,4,2,1]
	v[1,5] <- s[1,1,2,2] + s[4,4,3,3] - s[1,1,3,3] - s[4,4,2,2]
	v[2,3] <- s[1,2,1,3] + s[3,4,2,4] - s[1,2,2,4] - s[3,4,1,3]
	v[2,4] <- s[1,2,2,1] + s[3,4,4,3] - s[1,2,4,3] - s[3,4,2,1]
	v[2,5] <- s[1,2,2,2] + s[3,4,3,3] - s[1,2,3,3] - s[3,4,2,2]
	v[3,4] <- s[1,3,2,1] + s[2,4,4,3] - s[1,3,4,3] - s[2,4,2,1]
	v[3,5] <- s[1,3,2,2] + s[2,4,3,3] - s[1,3,3,3] - s[2,4,2,2]
	v[4,5] <- s[2,1,2,2] + s[4,3,3,3] - s[2,1,3,3] - s[4,3,2,2]
	v <- v + t(v)
	v[1,1] <- s[1,1,1,1] + s[4,4,4,4] - s[1,1,4,4] - s[4,4,1,1]
	v[2,2] <- s[1,2,1,2] + s[3,4,3,4] - s[1,2,3,4] - s[3,4,1,2]
	v[3,3] <- s[1,3,1,3] + s[2,4,2,4] - s[1,3,2,4] - s[2,4,1,3]
	v[4,4] <- s[2,1,2,1] + s[4,3,4,3] - s[2,1,4,3] - s[4,3,2,1]
	v[5,5] <- s[2,2,2,2] + s[3,3,3,3] - s[2,2,3,3] - s[3,3,2,2]
	list(v=v, cutoff=cutoff)
} #function
