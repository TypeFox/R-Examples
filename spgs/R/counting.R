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

#counting.R:  R functions for counting pairs, triples, quadruples and 
#cylinders of symbols in sequences, together with functions to convert 
#the count arrays/tables into vectors.

pair.counts <- function(x, case=c("lower", "upper", "as is"), circular=TRUE)
{
#check arguments
	if (!is.character(x))
		x <- as.character(x)
#		stop("x must be a character vector or an object of class SeqFastadna.")
	case <- match.arg(case)
	if (case=="lower") x <- tolower(x)
	if (case=="upper") x <- toupper(x)
	if (!is.logical(circular)) circular <- TRUE
#Convert sequence to a sequence of values, 1, ..., <numberOfSymbols
  uniqueSymbols <- as.character(unique(x))
  numSymbols <- length(uniqueSymbols)
  tr <- 1:numSymbols
  names(tr) <- sort(uniqueSymbols)
  x <- tr[x]
#Count consecutive pairs in the sequence
	counts <- .C("PairCounts",
		as.integer(x),
		as.integer(length(x)),
		as.integer(numSymbols),
		as.integer(circular),
		counts=integer(numSymbols*numSymbols),
		PACKAGE="spgs"
	)$counts
	dim(counts) <- c(numSymbols, numSymbols)
	dimnames(counts) <- list(names(tr), names(tr))
	counts #return array of counts
} #function

triple.counts <- function(x, case=c("lower", "upper", "as is"), circular=TRUE)
{
#check arguments
if (!is.character(x))
	x <- as.character(x)
#			stop("x must be a character vector or an object of class SeqFastadna.")
	case <- match.arg(case)
	if (case=="lower") x <- tolower(x)
	if (case=="upper") x <- toupper(x)
	if (!is.logical(circular)) circular <- TRUE
#Convert sequence to a sequence of values, 1, ..., <numberOfSymbols>
	uniqueSymbols <- as.character(unique(x))
	numSymbols <- length(uniqueSymbols)
	tr <- 1:numSymbols
	names(tr) <- sort(uniqueSymbols)
	x <- tr[x]
#Count consecutive triples in the sequence
	counts <- .C("TripleCounts",
		as.integer(x),
		as.integer(length(x)),
		as.integer(numSymbols),
		as.integer(circular),
		counts=integer(numSymbols*numSymbols*numSymbols),
		PACKAGE="spgs"
	)$counts
	dim(counts) <- c(numSymbols, numSymbols, numSymbols)
	dimnames(counts) <- list(names(tr), names(tr), names(tr))
	counts #return array of counts
} #function

quadruple.counts <- function(x, case=c("lower", "upper", "as is"), circular=TRUE)
{
#check arguments
if (!is.character(x))
	x <- as.character(x)
#			stop("x must be a character vector or an object of class SeqFastadna.")
	case <- match.arg(case)
	if (case=="lower") x <- tolower(x)
	if (case=="upper") x <- toupper(x)
	if (!is.logical(circular)) circular <- TRUE
#Convert sequence to a sequence of values, 1, ..., <numberOfSymbols>
	uniqueSymbols <- as.character(unique(x))
	numSymbols <- length(uniqueSymbols)
	tr <- 1:numSymbols
	names(tr) <- sort(uniqueSymbols)
	x <- tr[x]
#Count consecutive quadruples in the sequence
	counts <- .C("QuadrupleCounts",
		as.integer(x),
		as.integer(length(x)),
		as.integer(numSymbols),
		as.integer(circular),
		counts=integer(numSymbols*numSymbols*numSymbols*numSymbols),
		PACKAGE="spgs"
	)$counts
#	counts <- table(x, c(x[2:n], x[1]), c(x[3:n], x[1:2]), c(x[4:n], x[1:3]))
	dim(counts) <- c(numSymbols, numSymbols, numSymbols, numSymbols)
	dimnames(counts) <- list(names(tr), names(tr), names(tr), names(tr))
  counts #return array of counts
} #function

	cyl2lag2.counts <- function(x, last.lag, case=c("lower", "upper", "as is"), circular=TRUE)
{
#check arguments
	if (!is.character(x))
		x <- as.character(x)
#		stop("x must be a character vector or an object of class SeqFastadna.")
	if (!is.numeric(last.lag) || length(last.lag)>1 || floor(last.lag)!=last.lag || last.lag<0)
		stop("last.lag must be a positive integer value")
	case <- match.arg(case)
	if (case=="lower") x <- tolower(x)
	if (case=="upper") x <- toupper(x)
	if (!is.logical(circular)) circular <- TRUE
#Convert sequence to a sequence of values, 1, ..., <numberOfSymbols>
	uniqueSymbols <- as.character(unique(x))
	numSymbols <- length(uniqueSymbols)
	tr <- 1:numSymbols
	symbols <- sort(uniqueSymbols)
	names(tr) <- symbols
	x <- tr[x]
#Count consecutive quadruples in the sequence
	counts <- .C("Cyl2lag2Counts",
		as.integer(x),
		as.integer(length(x)),
		as.integer(last.lag),
		as.integer(circular),
		quads=integer(numSymbols*numSymbols*numSymbols*numSymbols*(last.lag+1)),
		PACKAGE="spgs"
	)$quads
	dim(counts) <- c(numSymbols, numSymbols, numSymbols, numSymbols, last.lag+1)
	dimnames(counts) <- list(symbols, symbols, symbols, symbols, 0:last.lag)
  counts #return array of counts
} #function


		cylinder.counts <- function(x, cylinder, case=c("lower", "upper", "as is"), 
circular=TRUE)
{
#Check arguments
	if (!is.numeric(cylinder) || !is.vector(cylinder))
		stop("cylinder must be a numeric vector.")
	if (any(cylinder<=0 | cylinder!=floor(cylinder)))
		stop("cylinder must be a vector of positive whole numbers")
	cylinder <- sort(as.integer(cylinder)) #make sure the elements of cylinder appear in increasing order
	if (any(diff(cylinder)==0))
	  stop("cylinder must not have repeated entries")
#	n <- length(x)
	cylinderLen <- length(cylinder)
	case <- match.arg(case)
	if (case=="lower") x <- tolower(x)
	if (case=="upper") x <- toupper(x)
	if (!is.logical(circular)) circular <- TRUE
#Convert sequence to a sequence of values, 1, ..., <numberOfSymbols>
	uniqueSymbols <- as.character(unique(x))
	numSymbols <- length(uniqueSymbols)
	tr <- 1:numSymbols
	names(tr) <- sort(uniqueSymbols)
	x <- tr[x]
#Count consecutive quadruples in the sequence
	counts <- .C("CylinderCounts",
		as.integer(x),
		as.integer(length(x)),
		as.integer(cylinder),
		as.integer(cylinderLen),
		as.integer(numSymbols),
		as.integer(circular),
		counts=integer(numSymbols**cylinderLen),
		PACKAGE="spgs"
	)$counts
	dim(counts) <- rep(numSymbols, cylinderLen)
	dimnames(counts) <- rep(list(names(tr)), cylinderLen)
#	names(dimnames(counts)) <- dnn
	counts #return array of counts
} #function


array2vector <- table2vector <- function(x, sep=".", sort=FALSE, 
		rev=FALSE, ...)
{
#check arguments
	if (is.null(x) || length(x)==0) return(as.vector(x))
	if (!is.character(sep) || length(sep)!=1)
		stop("sep must be a character string")
#Fill in missing dimension names with default values
	dimnames(x) <- lapply(1:length(dim(x)), function(i) {
		if (length(dimnames(x)[[i]])==0)
			formatC(1:dim(x)[i], width=nchar(as.character(dim(x)[i])), flag="0")
		else
			dimnames(x)[[i]]
	})
#Get array indices of all elements in x
	idx <- arrayInd(1:length(x), .dim=dim(x), .dimnames=dimnames(x))
#Convert the array indices to their corresponding array dimension names
	tags <- lapply(1:ncol(idx), function(j) dimnames(x)[[j]][idx[,j]])
#Make a list holding everything necessary to call paste to concatenate the dimension names and obtain a sensible name for each array element
	func.args <- c(list(paste), tags, list(sep=sep))
#	for (i in 1:length(tags))
#		func.args[[1+i]] <- tags[[i]]
#	func.args[["sep"]] <- sep
	vec <- as.vector(x) #Convert x to a vector, which strips off the array dimension names
	names(vec) <- eval(as.call(func.args))
#Sort the vector by the names attribute if sort is TRUE
	if (sort)
	{
		if (rev) keys <- strreverse(names(vec)) #if rev is TRUE, sort on the reversed element names
		else keys <- names(vec) #otherwise, sort on the plain element names
		vec <- vec[order(keys, ...)]
	} #if
	vec
} #function
		