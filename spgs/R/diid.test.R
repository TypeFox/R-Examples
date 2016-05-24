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

diid.test <- function(x, type=c("lb", "ks"), method="holm",	lag=20, ...)
{
#Check arguments
	dname <- deparse(substitute(x)) #record name of data
	if (is.list(x)) stop("x may not be a list")
	if (!is.character(x)) #is x a character vector?
		x <- as.vector(x, "character") #no, so coerse it to mode character
	types <- c("lb.test", "diffsign.test", "turningpoint.test", "ks.unif.test", "chisq.unif.test", "rank.test")
	typeNames <- match.arg(type, choices=types, several.ok=TRUE)
	if (length(typeNames)!=length(type))
		stop("an invalid test was specified in type.")
	method <- match.arg(method, choices=p.adjust.methods, several.ok=FALSE)
	if (length(method)!=1)
		stop("An invalid correction method was specified in method.")
	if ("lb.test" %in% typeNames)
	{
		if (!is.numeric(lag) || length(lag)!=1 || floor(lag)!=lag || lag<1)
			stop("lag must be a positive integer.")
	} #if

#Fit Bernoulli scheme and generate possible disturbance
	diid <- diid.disturbance(x, estimates=TRUE, ...)
  
 	#Carry out tests
	tests <- lapply(typeNames, function(test) {
		if (test=="lb.test") eval(call(test, diid$disturbance, lag=lag))
		else eval(call(test, diid$disturbance))
	})

#Pack return
	rval <- list()
	rval$method <- "Composite test for a Bernoulli scheme"
	rval$statistics <- sapply(tests, function(h) h$statistic)
	rval$parameters <- sapply(tests, function(h) ifelse(!is.null(h$parameter), h$parameter, NA))
	rval$p.values <- sapply(tests, function(h) h$p.value)
	rval$adjusted.p.values <- p.adjust(rval$p.values, method)
	rval$methods <- sapply(tests, function(h) h$method)
	rval$data.name <- dname
	rval$adjust.method <- method
	rval$estimate <- diid$stat.dist
	class(rval) <- c("multiplehtests")
	rval
} #function


diid.disturbance <- function(x, random=TRUE, estimates=FALSE)
{
#Check arguments
  if (!is.character(x)) #is x a character vector?
    x <- as.vector(x, "character") #no, so coerse it to mode character
  if (length(random)!=1 || (!is.numeric(random) && !is.logical(random)))
    stop("random must be a single logical or numeric value.")
  noise.value <- 0.5 #default value for fixed, constant noise source
  if (is.numeric(random))
  {
    if (random<0 || random>1)
      stop("random must be in the range 0-1.")
    noise.value <- random
    random <- FALSE
  }
  n <- length(x)
  if (n<2)
    stop("Not enough data.")
#Estimate stationary distribution for an i.i.d. sequence that may have produced x
  counts <- c(table(x)) #count symbols and convert table of counts to named vector
  pi <- counts/sum(counts)
#Initialise cdf of stationary distribution
  iid.cdf <- c(0, cumsum(pi[-length(pi)])) #precalculate cdf
  names(iid.cdf) <- names(pi) #correctly label the elements
#Generate noise
  if (random) noise <- runif(n) #Unif(0,1) noise
  else noise <- rep(noise.value, n) #constant noise
#Generate feasible disturbance from observed sample path and estimated stationary distribution
  disturbance <- iid.cdf[x] + pi[x]*noise
  names(disturbance) <- NULL #remove names
  if (estimates)
    return(list(disturbance=disturbance, stat.dist=pi))
  else
    return(disturbance)
  } #function
