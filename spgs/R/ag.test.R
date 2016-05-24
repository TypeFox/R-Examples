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

ag.test<- function(x, type=c("interval", "simplex"))
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
	type <- match.arg(type)

#Compute test statistic and perform test
	stat <- abs(probs[1]+probs[3]-0.5)
	names(stat) <- "etaV*"
#Compute p-value
	if (type=="simplex")
	{
		p.value <- pagtestonsimplex(stat)
		method <- "Test of relative purine/pyrimidine frequency equivalence based on relative purine frequency and the 3-simplex"
#		method <- "Test of Nucleotide Frequency Compliance with A+G=C+T on the Three-Simplex\n"
	}
	else
	{
		p.value <- pagtestoninterval(stat)
		method <- "Test of relative purine/pyrimidine frequency equivalence based on relative purine frequency over the unit interval"
#		method <- "Test of Nucleotide Frequency Compliance with A+G=C+T on the Unit Interval\n"
	} #if
	names(p.value) <- "p-value"

#Return result
	rval <- list(statistic=stat, p.value=p.value, method=method,
		data.name=dname, estimate=probs,
		stat.desc= "etaV* = abs(A+G-0.5)",
		null="A+G != C+T",
		alternative="A+G = C+T")
	class(rval) <- "htest.ext"
	rval
} #function

pagtestonsimplex <- function(stat)
#The cumulative distribution function of etaV* = |A+G-0.5| based on uniform 
#distribution over the 3-simplex.
{
	dist <- 3*stat - 4*stat*stat*stat
	dist[stat<0] <- 0
	dist[stat>0.5] <- 1
	dist
} #function

pagtestoninterval <- function(stat)
#The cumulative distribution function of etaV* = |A+G-0.5| based on uniform 
#distribution over the unit interval.
{
	dist <- 2*stat
	dist[stat<0] <- 0
	dist[stat>0.5] <- 1
	dist
} #function
