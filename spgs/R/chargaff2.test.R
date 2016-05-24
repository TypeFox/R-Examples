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

chargaff2.test <- function(x, alg=c("table", "simulate", "upper"), n, no.p.value=FALSE)
{
#Check arguments
	dname <- deparse(substitute(x))
	if (is.character(x) || class(x)=="SeqFastadna")
	{ #extract relative frequencies of nucleic acids from DNA sequence
		probs <- pair.counts(x)
		x<- probs/rowSums(probs) #normalise row sums to unity
		x[is.na(x)] <- 0.0 #fix rows of zeros
	}
#	else
#		x <- as.matrix(x) #convert to matrix, leaving names in tact
	probs <- ISPX2mat(x) #prepare correctly ordered 4X4 stochastic matrix
	alg <- match.arg(alg)
	if (alg=="simulate")
	{
		if (missing(n)) stop("n has not been specified")
		if (!is.numeric(n) || length(n)!=1 || n<=0 || n!=floor(n)) 
  		stop("n must be a positive integer")
  } #if

#Compute test statistic and perform test
	probs24 <- 1-(probs[2,1]+probs[2,2]+probs[2,3])
	f1 <- probs[2,2]
	f2 <- 1 - probs[3,1] - probs[2,2] - probs24*probs[1,2]/probs[1,3]
	f3 <- probs[2,1]*probs[1,3]/probs24
	f4 <- probs[3,1]*probs[1,3]/probs24
	f5 <- 1 - probs[1,1] - (probs[2,1]+probs[3,1])*probs[1,3]/probs24
	raw = c(f1-probs[3,3], f2-probs[3,2], f3-probs[4,3], f4-probs[4,2], f5-probs[4,1])
	if (f2>0 && f5>0) stat <- max(abs(raw))
	else stat <- 1
	names(stat) <- "eta2"

#if no.p.value is FALSE, compute p-value using the specified method
	if (!no.p.value)
	{
		p.value <- switch(alg,
			table=chargaff2testPValueFromTable(stat),
			simulate=chargaff2testPValueFromMonteCarlo(stat, n),
			upper=chargaff2testPValueFromUpperBound(stat)
		) #switch
		names(p.value) <- "p-value"
	} #if
	
#Return result
	method <- "Matrix test of Chargaff's second parity rule (CSPR) for dinucleotides\n"
	method <- switch(alg, 
		table=paste(method, "with p-value from linearly interpolated look-up table"),
		simulate=paste(method, "with simulated p-value (based on ",n, " replicates)"),
		lower=paste(method, "with lower bound computed for p-value"),
		upper=paste(method, "with upper bound computed for p-value")
	) #switch
	rval <- list(statistic=stat, method=method,
		data.name=dname, f=c(f1=f1, f2=f2, f3=f3, f4=f4, f5=f5), estimate=probs,
		null="The DNA sequence (or matrix P) does not comply with CSPR for dinucleotides",
		alternative="The DNA sequence (or matrix P) complies with CSPR for dinucleotides")
	if (!no.p.value) rval$p.value <- p.value
	class(rval) <- "htest.ext"
	rval
} #function

chargaff2testPValueFromTable <- function(stat)
{
	if (!is.numeric(stat)) stop("stat must be numeric")
	approx(ch2$quantile, ch2$prob, as.vector(stat), method="linear", rule=2)$y
} #function

chargaff2testPValueFromMonteCarlo <- function(epsilon, n)
{
	if (missing(epsilon)) stop("epsilon has not been specified")
	if (!is.numeric(epsilon) || length(epsilon)!=1) stop("epsilon must be a single numeric value")
	if (epsilon<=0) return(0)
	if (epsilon>=1) return(1)
	if (missing(n)) stop("n has not been specified")
	if (!is.numeric(n) || length(n)!=1 || n<=0 || n!=floor(n)) 
  	stop("n must be a positive integer")
	p <- rexp(16*n, 1)
	p.val <- .C("ComputeEta2Statistic", as.double(p), as.integer(n), as.double(epsilon),
		res=double(1), PACKAGE="spgs")$res
	p.val
} #function

chargaff2testPValueFromMonteCarlo2 <- function(epsilon, n, diag=FALSE)
{
	if (missing(epsilon)) stop("epsilon has not been specified")
	if (!is.numeric(epsilon) || length(epsilon)!=1) stop("epsilon must be a single numeric value")
	if (epsilon<=0) return(0)
	if (epsilon>=1) return(1)
	if (missing(n)) stop("n has not been specified")
	if (!is.numeric(n) || length(n)!=1 || n<=0 || n!=floor(n)) 
  	stop("n must be a positive integer")
#p <- array(rexp(16*n, 1), c(n, 4, 4))
#x <- p/array(rep(apply(p, c(1,2), sum), 4), dim=dim(p))
	p <- rexp(16*n, 1)
	x <- .C("ProbabilityNormalise", as.double(p), as.integer(n), as.integer(4), 
		as.integer(4), res=double(16*n),
		PACKAGE="spgs")$res
	x <- array(x, dim=c(n, 4, 4))
	x24 <- 1-(x[,2,1]+x[,2,2]+x[,2,3])
	f1 <- x[,2,2]
	f2 <- 1 - x[,3,1] - x[,2,2] - x24*x[,1,2]/x[,1,3]
	f3 <- x[,2,1]*x[,1,3]/x24
	f4 <- x[,3,1]*x[,1,3]/x24
	f5 <- 1 - x[,1,1] - (x[,2,1]+x[,3,1])*x[,1,3]/x24
	raw = cbind(f1-x[,3,3], f2-x[,3,2], f3-x[,4,3], f4-x[,4,2], f5-x[,4,1])
	stat <- rep(1, n)
	good <- f2>0 & f5>0
	stat[good] <- apply(abs(raw[good,]), 1, max)
 p.val <- sum(stat<=epsilon)/n
	if (diag) list(p=p, x=x, x24=x24, f1=f1, f2=f2, f3=f3, f4=f4, f5=f5, raw=raw, good=good, stat=stat, p.value=p.val)
	else p.val
} #function

chargaff2testPValueFromUpperBound <- function(epsilon)
{
	min(c(27*epsilon*epsilon, 1))
} #function
