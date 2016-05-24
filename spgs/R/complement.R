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

#complement.R

complement <- function(x, ...)
	UseMethod("complement")

complement.default <- function(x, content=c("dna", "rna"), 
case=c("lower", "upper", "as is"), ...)
#returns the complement of the DNA or RNA sequences in x.
{
#Check arguments
	if (!is.character(x)) x <- as.character(x) #coerse x to a character vector
	content <- match.arg(content)
	case <- match.arg(case)
	if (length(x)==0 || (length(x)==1 && nchar(x)==0))
		return(x) #bail if input is empty
	if (case=="lower") x <- tolower(x)
	if (case=="upper") x <- toupper(x)
	if (content=="dna")
	{
		src  <- "acgturykmswbdhvnxACGTURYKMSWBDHVNX-"
		dest <- "tgcaayrmkswvhdbnxTGCAAyRMKSWVHDBNX-"
	}
	else
	{
		src  <- "acgturykmswbdhvnxACGTURYKMSWBDHVNX-"
		dest <- "ugcaayrmkswvhdbnxUGCAAyRMKSWVHDBNX-"
	} #if
	chartr(src, dest, x)
} #function

complement.SeqFastadna <- function(x, ...)
{
	seq <- complement.default(as.character(x), ...)
	attributes(seq) <- attributes(x)
	seq
} #function

complement.list <- function(x, ...)
{
	lapply(x, function(y) complement(y, ...))
} #function
