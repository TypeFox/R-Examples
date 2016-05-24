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

FillToSquare <- function(x)
{
	if (is.null(rownames(x)) || is.null(colnames(x)))
		stop("x must have both column and row names")
	allNames <- sort(union(rownames(x), colnames(x)))
	missing.rows <- setdiff(allNames, rownames(x))
	missing.cols <- setdiff(allNames, colnames(x))
	x <- rbind(x, matrix(0, nrow=length(missing.rows), ncol=ncol(x), dimnames=list(missing.rows, colnames(x))))
	x <- cbind(x, matrix(0, nrow=nrow(x), ncol=length(missing.cols), dimnames=list(rownames(x), missing.cols)))
	x[allNames, allNames]
} #function

statdist <- function(x, copy.names=TRUE)
{
	if (is.list(x) || class(x)=="FiniteStateMarkovChain")
	{
		if (is.null(x$stat.dist))
			stop("x is not a valid list or FiniteStateMarkovChain class")
		else
			return(x$stat.dist)
	} #if
	if (!is.matrix(x) || nrow(x)!=ncol(x)) stop("x must be a square matrix.")
	if (any(x<0) ||!all.equal(rowSums(x), rep(1, nrow(x)), check.attributes=FALSE))
		stop("x must contain non-negative elements and have row sums equal to 1.")
	if (copy.names && !identical(rownames(x), colnames(x)))
		stop("The row and column names of the transition matrix are not consistent.")
	e <- eigen(t(x))
	stat.dist <- abs(as.double(e$vectors[, which.max(abs(e$values))]))
	stat.dist <- stat.dist/sum(stat.dist)
	if (copy.names) names(stat.dist) <- colnames(x)
	stat.dist
} #function

estimateMarkovChain <- function(x, circular=TRUE)
{
	if (length(x)==0) stop("x is an empty sequence")
	if (!is.character(x)) x <- as.character(x)
#Estimate transition matrix and stationary distribution for a Markov chain that may have produced x
	states <- sort(unique(x))
	if (circular)
  	counts <- table(x, c(x[-1], x[1]), dnn=NULL)
  else
  	counts <- table(x[-length(x)], x[-1], dnn=NULL)
  counts <- FillToSquare(counts) #insert 0 counts for transitions that never occur in x
  pi <- c(rowSums(counts))
  p <- as.matrix(counts/pi)
  pi <- pi/sum(pi)
  p[is.na(p)] <- 0
  chain <- list(trans.mat=p, stat.dist=pi, states=states)
  class(chain) <- "FiniteStateMarkovChain"
  chain
} #function

print.FiniteStateMarkovChain <- function(x, ...)
{
	cat("Finite State Markov Chain\nNumber of states:", length(x$states), "\n")
	cat("States:", x$states, "\n")
	cat("Stationary distribution:\n")
	print(x$stat.dist)
	cat("Transition matrix:\n")
	print(x$trans.mat)
} #function

simulateMarkovChain <- function(n, trans.mat, init.dist=NULL, states=colnames(trans.mat))
{
#Check arguments
	if (!is.numeric(n) || !is.vector(n) || length(n)!=1 || floor(n)!=n || n<1)
	stop("n must be a positive integer.")
	if (!is.array(trans.mat) || length(dim(trans.mat))!=2 || nrow(trans.mat)!=ncol(trans.mat))
		stop("trans.mat must be a square, two dimensional array, for example, a matrix or a table")
	if (is.null(states))
		if (!is.null(rownames(trans.mat))) states <- rownames(trans.mat)
		else states <- 1:dim(trans.mat)[2]
	if (is.null(init.dist))
	{
		e <- eigen(t(trans.mat))
		init.dist <- abs(as.double(e$vectors[,1]))
		init.dist <- init.dist/sum(init.dist)
	} #if
	if (length(init.dist)!=ncol(trans.mat))
		stop("The dimensions of trans.mat and length of init.dist do not match")
	disturbance <- runif(n) #generate U(0,1) noise
	samplePath <- .C("GenerateMarkovSamplePath", #function name
		as.double(trans.mat), #transition matrix of Markov chain
		as.double(init.dist), #probability function of initial state
		as.integer(length(states)), #number of states
		as.double(disturbance), #source of U(0,1) random noise to drive the simulation
		as.integer(n), #length of sample path to simulate
		path=integer(n), #workspace to hold simulation
		PACKAGE="spgs"
	)$path
	states[samplePath] #convert numeric simulation to specified state space representation and return
} #function
	