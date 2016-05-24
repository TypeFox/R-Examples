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

agct.test <- function(x, alg=c("exact", "simulate", "lower", "Lower", "upper"), 
n)
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
  		stop("n must be a positive integer")
  } #if

#Compute test statistic and perform test
	xstar <- 0.25+0.5*(probs[1]-probs[3])
	if (xstar<0) xstar <- 0
	if (xstar>0.5) xstar <- 0.5
	ystar <- 0.25+0.5*(probs[2]-probs[4])
	if (ystar<0) ystar <- 0
	if (ystar>0.5) ystar <- 0.5
	offset <- probs-c(xstar, ystar, 0.5-xstar, 0.5-ystar)
	stat <- c(etaV=sqrt(sum(offset*offset)))

#Compute p-value using the specified method
	p.value <- switch(alg,
		exact=agctTestPValue(stat),
		simulate=agctTestPValueFromMonteCarlo(stat, n),
		lower=agctTestPValueLowerBound(stat),
		Lower=agctTestPValueLowerBoundTight(stat),
		upper=agctTestPValueUpperBound(stat)
	) #switch
	names(p.value) <- "p-value"

 #Return result
	method <- "Test of relative purine/pyrimidine frequency equivalence based on Euclidean distance\n"
	method2 <- switch(alg, 
		exact="",
		simulate=sprintf(" with simulated p-value (based on %d replicates)", n),
		lower=" with lower bound computed for p-value",
		Lower=" with alternative, tighter lower bound computed for p-value",
		upper=" with upper bound computed for p-value"
	) #switch
	method <- sprintf("%s%s", method, method2)
	rval <- list(statistic=stat, p.value=p.value, method=method,
		data.name=dname, estimate=probs, 
		stat.desc= "etaV = Euclidean distance from relative frequency vector to closest point in\nthetaV ={(x,y,1/2-x,1/2-y) : 0 <= x,y <= 1/2}",
		null="A+G != C+T",
		alternative="A+G = C+T")
	class(rval) <- "htest.ext"
	rval
} #function

pagcttest <- function(stat)
#The cumulative distribution function of etaV based on uniform distribution 
#over the 3-simplex.
{
	if (missing(stat)) stop("stat has not been specified")
	if (!is.vector(stat, "numeric") || length(stat)==0)
		stop("stat must be a non-empty numeric vector")
	sapply(stat, agctTestPValue)
} #function

agctTestPValue <- function(epsilon)
{
	if (epsilon<0) return(0)
	if (epsilon>sqrt(3/8)) return(1)
#Declare function for computing area of trapezoids in the corners of the simplex
	intCurvedCornerTrapezoids <- function(epsilon, x)
	{
		epsilon*(epsilon*(asin(x)+x*sqrt(1-x*x))/2 - epsilon*x*x/sqrt(2) + x*x*x*epsilon*epsilon-epsilon*epsilon*x)/sqrt(2)
	} #function

	if (epsilon<=0.5)
	{
		triPrismVol <- 3*epsilon*(1-epsilon)
		shrinkingTriangleVol <- 0
		curvedVol <- 12*(intCurvedCornerTrapezoids(epsilon, 1/sqrt(3)) - intCurvedCornerTrapezoids(epsilon, 0))
	}
	else
	{
	 triPrismVol <- 0.75
	 shrinkingTriangleVol <- 0.25 - 2*(0.5-sqrt(2*epsilon*epsilon-.5))^3
	 curvedVol <- 12*(intCurvedCornerTrapezoids(epsilon, 1/sqrt(3)) - intCurvedCornerTrapezoids(epsilon, sqrt(1-0.25/(epsilon*epsilon))))
	} #if
	triPrismVol+shrinkingTriangleVol +curvedVol
} #function

agctTestPValueFromMonteCarlo <- function(epsilon, n)
{
	if (missing(epsilon)) stop("epsilon has not been specified")
	if (!is.numeric(epsilon) || length(epsilon)!=1) stop("epsilon must be a single numerical value")
	if (epsilon<=0) return(0)
	if (epsilon>=1) return(1)
	if (missing(n)) stop("n has not been specified")
	if (!is.numeric(n) || length(n)!=1 || n<=0 || n!=floor(n)) 
  	stop("n must be a positive integer")
	rexp <- rexp(4*n, 1)
	pr <- .C("ProbabilityNormalise", as.double(rexp), as.integer(n), as.integer(1), 
		as.integer(4), res=double(4*n), PACKAGE="spgs")$res
	dim(pr) <- c(n,4)
	xstar <- 0.25+0.5*(pr[,1]-pr[,3])
	xstar[xstar<0] <- 0
	xstar[xstar>0.5] <- 0.5
	ystar <- 0.25+0.5*(pr[,2]-pr[,4])
	ystar[ystar<0] <- 0
	ystar[ystar>0.5] <- 0.5
	offset <- pr-cbind(xstar, ystar, 0.5-xstar, 0.5-ystar)
	stat <- sqrt(rowSums(offset*offset))
	return(sum(stat<=epsilon)/n)
} #function

agctTestPValueLowerBound <- function(epsilon)
{
	3*epsilon*(1-epsilon)
} #function

agctTestPValueLowerBoundTight <- function(epsilon)
{
	l <- sqrt(2/3)*epsilon
#	3*l-4*l^3 + 3*(epsilon-l)*(1-l-epsilon)
	3*l-4*l^3 + 3*(epsilon*(1-epsilon)-l*(1-l))
} #function

agctTestPValueUpperBound <- function(epsilon)
{
	3*epsilon-4*epsilon^3
} #function
