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

#disambiguate.R

charfilter <- function(x, symbols)
#Eliminates all but the specified symbols from character strings in a vector.
{
	if (!is.character(x))
		stop("x must be a character vector")
	if (length(symbols)>1 && max(nchar(symbols))>1)
		stop("symbols must be a single character string or a vector of single character elements")
	if (length(symbols)==1 && nchar(symbols)>1)
		symbols <- split(symbols, "")
	sapply(strsplit(x, ""), function(y) paste(y[y %in% symbols], collapse=""))
} #function

disambiguate <- function(x, ...)
	UseMethod("disambiguate")

disambiguate.default <- function(x, case=c("lower", "upper", "as is"), ...)
{
#Check arguments
	if (!is.character(x)) x <- as.character(x) #coerse x to a character vector
	case <- match.arg(case)
	if (length(x)==0 || (length(x)==1 && nchar(x)==0))
		return(x) #bail if input is empty
	if (case=="lower")
	{
		x <- tolower(x)
		symbols <- c("a", "c", "g", "t", "u")
	}
	else if (case=="upper")
	{
		x <- toupper(x)
		symbols <- c("A", "C", "G", "T", "U")
	}
	else
		symbols <- c("a", "c", "g", "t", "u", "A", "C", "G", "T", "U")
	if (max(nchar(x))>1) #is x a vector of one or more non-single character strings?
		return(charfilter(x, symbols))
#x is not a single string, so process it as a vector
	x[x %in% symbols]
} #function

disambiguate.SeqFastadna <- function(x, ...)
{
	seq <- disambiguate.default(as.character(x), ...)
	attributes(seq) <- attributes(x)
	seq
} #function

disambiguate.list <- function(x, ...)
{
	lapply(x, function(y) disambiguate(y, ...))
} #function
