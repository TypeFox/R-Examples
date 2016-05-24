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

rstochmat <- function(n, labels)
{
	if (missing(n))
	{
		if (missing(labels) || length(labels)==0)
			stop("No valid arguments have been specified")
		else
			n <- length(labels)
	} #if
	if (!is.numeric(n) || length(n)!=1 || is.infinite(n) || n<1 || floor(n)!=n)
		stop("n must be a positive integer value")
	if (missing(labels)) labels <- as.character(1:n)
  mat <- matrix(rexp(n*n), n, n, dimnames=list(labels, labels))
  mat/rowSums(mat) #make the matrix stochastic
} #function

rcspr2mat <- function(labels=c("a", "c", "g", "t"))
{
n <- matrix(0, 4, 4)
n[1,1] <- -log(runif(1))
n[2,2] <- -log(runif(1))
n[3,3] <- n[2,2]
n[4,4] <- n[1,1]

gamma3i <- 0
s <- Inf
t <- Inf
while (gamma3i<=s || gamma3i<=t)
{
  gamma3o <- -sum(log(runif(3)))
  gamma3i <- -sum(log(runif(3)))

uu <- gamma3o*runif(2)
u <- diff(c(0, min(uu), max(uu), gamma3o))
vv <- gamma3o*runif(2)
v <- diff(c(0, min(vv), max(vv), gamma3o))
n[2:4, 1] <- u
n[1, 2:4] <- v
n[2:3, 4] <- v[2:1]
n[4, 2:3] <- u[2:1]

  s <- u[1]+v[2]
  t <- u[2]+v[1]
#while gamma3i<=s || gamma3i<=t
#  gamma3i <- -sum(log(rand(3, 1)));#
#end %while

  if (gamma3i>s && gamma3i>t)
  {
    n[2,3] <- gamma3i-s
    n[3,2] <- gamma3i-t
  } #if
} #while

if (runif(1)<=0.5) n <- t(n)
if (runif(1)<=0.5) n <- n[c(1, 3, 2, 4), c(1, 3, 2, 4)]
if (runif(1)<=0.5) n <- n[c(4, 2, 3, 1), c(4, 2, 3, 1)]

dimnames(n) <- list(labels, labels)
n/rowSums(n)
} #function

rstochvec <- function(n, labels)
{
	if (missing(n))
	{
		if (missing(labels) || length(labels)==0)
			stop("No valid arguments have been specified")
		else
			n <- length(labels)
	} #if
	if (!is.numeric(n) || length(n)!=1 || is.infinite(n) || n<1 || floor(n)!=n)
		stop("n must be a positive integer value")
	if (missing(labels)) labels <- as.character(1:n)
  vec <- rexp(n)
  names(vec) <- labels
  vec/sum(vec) #make the vector stochastic
} #function

