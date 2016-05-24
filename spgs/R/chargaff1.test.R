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

chargaff1.test <- function(x, alg=c("table", "simulate", "upper"), n, no.p.value=FALSE)
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
	f1 <- (1-(probs[1,1]+probs[4,1]))*(1-(probs[2,2]+probs[3,2]))/(probs[2,1]+probs[3,1])-probs[1,2]
	f2 <- (1-(probs[2,3]+probs[3,3]))/(1-(probs[2,2]+probs[3,2])) * (probs[1,2]+f1) - probs[1,3]
	raw <- c(f1-probs[4,2], f2-probs[4,		3])
	h <- (probs[1,1]+probs[4,1]<1) && (probs[2,2]+probs[3,2]<1) && (probs[2,3]+probs[3,3]<1)
	h <- h && (f1>0) && (f2>0) && (probs[4,1]+f1+f2<1)
	if (h) stat <- max(abs(raw))
	else stat <- 1
	names(stat) <- "eta1"

#if no.p.value is FALSE, compute p-value using the specified method
	if (!no.p.value)
	{
		p.value <- switch(alg,
			table=chargaff1testPValueFromTable(stat),
			simulate=chargaff1testPValueFromMonteCarlo(stat, n),
			upper=chargaff1testPValueFromUpperBound(stat)
		) #switch
		names(p.value) <- "p-value"
	} #if

#Return result
	method <- "Matrix test of Chargaff's second parity rule (CSPR) for mononucleotides\n"
	method <- switch(alg, 
		table=paste(method, "with p-value from linearly interpolated look-up table"),
		simulate=paste(method, "with simulated p-value (based on ",n, " replicates)"),
		lower=paste(method, "with lower bound computed for p-value"),
		upper=paste(method, "with upper bound computed for p-value")
	) #switch
	rval <- list(statistic=stat, method=method,
		data.name=dname, f=c(f1=f1, f2=f2), estimate=probs,
		null="A != T or C != G",
		alternative="A = T and C = G")
	if (!no.p.value) rval$p.value <- p.value
	class(rval) <- "htest.ext"
	rval
} #function

chargaff1testPValueFromTable <- function(stat)
{
	if (!is.numeric(stat)) stop("stat must be numeric")
	approx(ch1$quantile, ch1$prob, as.vector(stat), method="linear", rule=2)$y
} #function

chargaff1testPValueFromMonteCarlo <- function(epsilon, n)
{
	if (missing(epsilon)) stop("epsilon has not been specified")
	if (!is.numeric(epsilon) || length(epsilon)!=1) stop("epsilon must be a single numeric value")
	if (epsilon<=0) return(0)
	if (epsilon>=1) return(1)
	if (missing(n)) stop("n has not been specified")
	if (!is.numeric(n) || length(n)!=1 || n<=0 || n!=floor(n)) 
  	stop("n must be a positive integer")
	p <- rexp(16*n, 1)
	p.val <- .C("ComputeEta1Statistic", as.double(p), as.integer(n), as.double(epsilon),
		res=double(1), PACKAGE="spgs")$res
	p.val
} #function

chargaff1testPValueFromUpperBound <- function(epsilon)
{
	min(c(20*epsilon*epsilon/81, 1))
} #function
