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

markov.test <- function(x, type=c("lb", "ks"), method="holm",	lag=20, ...)
{
#Check arguments
	dname <- deparse(substitute(x))
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

#Estimate Markov chain and generate possible disturbance
	chain <- estimateMarkovChain(x) #estimate transition matrix and stationary distribution for x
		disturbance <- markov.disturbance(x, chain,estimates=FALSE, ...) #obtain a possible disturbance

#Carry out tests
	tests <- lapply(typeNames, function(test) {
		if (test=="lb.test") eval(call(test, disturbance, lag=lag))
		else eval(call(test, disturbance))
	})

#Pack return
	rval <- list()
	rval$method <- "Composite test for a first-order (finite state) Markov chain"
	rval$statistics <- sapply(tests, function(h) h$statistic)
	rval$parameters <- sapply(tests, function(h) ifelse(!is.null(h$parameter), h$parameter, NA))
	rval$p.values <- sapply(tests, function(h) h$p.value)
	rval$adjusted.p.values <- p.adjust(rval$p.values, method)
	rval$methods <- sapply(tests, function(h) h$method)
	rval$data.name <- dname
	rval$adjust.method <- method
	rval$estimate <- chain$trans.mat
	class(rval) <- c("multiplehtests")
	rval
} #function


markov.disturbance <- function(x, chain=NULL, random=TRUE, bandwidth=1, 
    estimates=is.null(chain))
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
  if (!is.numeric(bandwidth) || !is.vector(bandwidth) || length(bandwidth)!=1 || bandwidth<0 || bandwidth>1)
    stop("bandwidth must be a value in the range 0-1.")
  n <- length(x)
  if (n<2)
    stop("Not enough data.")
#If necessary, estimate transition matrix and stationary distribution for a Markov chain that may have produced x
  uniqueStates <- unique(x)
  if (is.null(chain))
    chain <- estimateMarkovChain(x)
  if (!("trans.mat" %in% names(chain))) 
  	stop("chain is lacking a transition matrix specification.")
  if (!("stat.dist" %in% names(chain)))
  	stop("chain is lacking a stationary distribution specification.")
  p <- chain$trans.mat
  pi <- chain$stat.dist
#Check state labels of p and pi
  if (!identical(rownames(p), colnames(p))) stop("The row and column names of the transition matrix are not consistent.")
  if (!identical(names(pi), colnames(p))) stop("The row and column names of the transition matrix are not consistent with the state labels on the stationary distribution.")
  if (is.null(names(pi))) stop("the stationary distribution is not labelled by states.")
  if (is.null(colnames(p))) stop("The rows and columns of the transition matrix are not labelled by states.")
#Initialise cdf's of initial state and conditional state transitions
  if (dim(p)[2L]==1)
  {
    disturbance <- rep(names(pi)[1L], n)
    if (estimates) list(disturbance=disturbance, trans.mat=chain$trans.mat, stat.dist=chain$stat.dist)
    else disturbance
  } #if
  init.cdf <- c(0, cumsum(pi[-length(pi)])) #precalculate cdf of initial state
  names(init.cdf) <- names(pi) #correctly label the elements
  if (dim(p)[2L]==2L) #if p is 2X2
    trans.cdf <- matrix(c(0, 0, p[,2L]), 2L, 2L)
  else #p is larger than 2X2
    trans.cdf <- cbind(rep(0, dim(p)[1L]), t(apply(p[, -dim(p)[2L]], 1L, cumsum))) #precalculate state-dependent map from closed unit interval  to state space
  dimnames(trans.cdf) <- dimnames(p) #correctly label the rows and columns
#Generate noise
  if (random)
    noise <- 0.5+bandwidth*(runif(n)-0.5) #uniform jitter centred around 0.5
  else
    noise <- rep(noise.value, n) #constant noise
#Generate feasible disturbance from observed sample path and estimated transition matrix and stationary distribution
  disturbance1 <- init.cdf[x[1]] + pi[x[1]]*noise[1] #get realisation of noise giving rise to the initial state
  xInd <- cbind(x[-n], x[-1]) #build index array
  disturbance <- trans.cdf[xInd] + p[xInd]*noise[2:n] #get the disturbance giving rise to the i-th state
  disturbance <- c(disturbance1, disturbance)
  if (estimates) list(disturbance=disturbance, trans.mat=chain$trans.mat, stat.dist=chain$stat.dist)
  else disturbance
} #function
