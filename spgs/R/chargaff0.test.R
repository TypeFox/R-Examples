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

	chargaff0.test <- function(x, alg=c("exact", "simulate", "lower", "upper", "Lower", "Upper"),
n, no.p.value=FALSE)
{
#Check arguments
	dname <- deparse(substitute(x))
	if (is.character(x) || class(x)=="SeqFastadna")
	{ #extract relative frequencies of nucleic acids from DNA sequence
		probs <- table(x)
		x <- c(probs/sum(probs)) #probability-normalise vector
	}
	else
		x <- c(x) #convert to vector, leaving names in tact
	probs <- ISPX2vec(x) #prepare correctly ordered probability vector
	alg <- match.arg(alg)
	if (alg=="simulate")
	{
		if (missing(n)) stop("n has not been specified")
		if (!is.numeric(n) || length(n)!=1 || n<=0 || n!=floor(n)) 
  		stop("n must be a positive integer ")
  } #if

#Compute statistic and perform test
	papt <- probs[1]-probs[4]
	pcpg <- probs[2]-probs[3]
	stat <- sqrt(0.5*(papt*papt+pcpg*pcpg))
	names(stat) <- "eta0"

#if no.p.value is FALSE, compute p-value using the specified method
	if (!no.p.value)
	{
		p.value <- switch(alg,
			exact=chargaff0testPValue(stat),
			simulate=chargaff0testPValueFromMonteCarlo(stat, n),
			lower=chargaff0testPValueLowerBound(stat),
			upper=chargaff0testPValueUpperBound(stat),
			Lower=chargaff0testPValueLowerBoundTight(stat),
			Upper=chargaff0testPValueUpperBoundTight(stat)
		) #switch
		names(p.value) <- "p-value"
	} #if

#Return result
	method <- "Vector test of Chargaff's second parity rule (CSPR) for mononucleotides"
	method <- switch(alg, 
		exact=method,
		simulate=paste(method, "with simulated p-value (based on ",n, " replicates)"),
		lower=paste(method, "with lower bound computed for p-value"),
		upper=paste(method, "with upper bound computed for p-value"),
		Lower=paste(method, "with alternative, tighter lower bound computed for p-value"),
		Upper=paste(method, "with alternative, tighter upper bound computed for p-value")
	) #switch
	rval <- list(statistic=stat, method=method,
		data.name=dname, estimate=probs,
		stat.desc= "eta0 = sqrt((A-T)^2/2 + (C-G)^2/2)",
		null="A != T or C != G",
		alternative="A = T and C = G")
	if (!no.p.value) rval$p.value <- p.value
	class(rval) <- "htest.ext"
	rval
} #function

chargaff0testPValueFromMonteCarlo <- function(epsilon, n)
{
	if (missing(epsilon)) stop("epsilon has not been specified")
	if (!is.numeric(epsilon) || length(epsilon)!=1) stop("epsilon must be a single numeric value")
	if (epsilon<=0) return(0)
	if (epsilon>0.25*sqrt(2)) return(1)
	if (missing(n)) stop("n has not been specified")
	if (!is.numeric(n) || length(n)!=1 || n<=0 || n!=floor(n)) 
  	stop("n must be a positive integer")
  pr <- matrix(runif(3*n), ncol=3)
  pr <- cbind(2*pr[, 1:2]-1, pr[, 3])
  cut <- sum(pr[,1]*pr[,1]+pr[,2]*pr[,2]<=1 & pr[,1]>pr[,3])
  cut <- cut/n*4*sqrt(2)*epsilon^3
  return(3*pi*epsilon^2 - 12*cut)
} #function

chargaff0testPValue <- function(epsilon)
{
	if (epsilon<0) return(rep(0, length(epsilon)))
	if (epsilon>0.25*sqrt(2)) return(rep(1, length(epsilon)))
	cylinderVolume <- 3*pi*epsilon^2*(1-2*sqrt(2)*epsilon)
	endVolume <- 2*sqrt(2)*(3*pi-4)*epsilon^3
	cylinderVolume+endVolume
} #function

#chargaff0testPValueFromTable <- function(stat)
#{
#	if (!is.numeric(stat)) stop("stat must be numeric")
#	pval <- approx(ch0$quantile, ch0$prob, as.vector(stat), method="linear", rule=2)$y
#	if (pval<ch0[2,1]) pval <- chargaff0testPValueUpperBoundTight(stat)
#	pval
#} #function

chargaff0testPValueLowerBound <- function(epsilon)
{
	if (epsilon<0) return(rep(0, length(epsilon)))
	if (epsilon>0.25*sqrt(2)) return(rep(1, length(epsilon)))
	3*pi*epsilon*epsilon*(1-2*sqrt(2)*epsilon)
} #function

chargaff0testPValueUpperBound <- function(epsilon)
{
		if (epsilon<0) return(rep(0, length(epsilon)))
	if (epsilon>0.25*sqrt(2)) return(rep(1, length(epsilon)))
3*pi*epsilon*epsilon
} #function

chargaff0testPValueLowerBoundTight <- function(epsilon)
{
	if (epsilon<0) return(rep(0, length(epsilon)))
	if (epsilon>0.25*sqrt(2)) return(rep(1, length(epsilon)))
	3*pi*epsilon*epsilon*(1-2*sqrt(2)*epsilon) + sqrt(2)*pi*epsilon*epsilon*epsilon
} #function

chargaff0testPValueUpperBoundTight <- function(epsilon)
{
	if (epsilon<0) return(rep(0, length(epsilon)))
	if (epsilon>0.25*sqrt(2)) return(rep(1, length(epsilon)))
	3*pi*epsilon*epsilon*(1-2*sqrt(2)*epsilon) + 12*sqrt(2)*epsilon*epsilon*epsilon
} #function
