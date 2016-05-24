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

diffsign.test <- function(x)
{
#check arguments
  if (!is.numeric(x)) #is x numeric?
    stop("x must be numeric.")
  if (NCOL(x)>1)
    stop("x is not a vector or univariate time series")
  dname <- deparse(substitute(x)) #name of data variable
  x <- c(x) #coerse to a vector and strip names
#Trim repeated values from the data
  x <- x[diff(x)!=0]
#Perform the test
  n <- length(x) #length of data set
  if (n<2)
    stop("there is not enough data to perform the difference-sign test")
  npos <- sum(diff(x)>0) #count positive differences
  mu <- (n-1)/2 #the expected number of positive differences in n IID observations
  sigma <- sqrt((n+1)/12) #the variance of the number of positive differences in n IID observations
  stat <- (npos-mu)/sigma #compute the test statistic which is asymptotically standard normal
  names(stat) <- "D" #name of statistic
  p <- 2*pnorm(abs(stat), lower.tail=FALSE) #get the p-value of the test statistic
#Pack return
  rval <- list(statistic=stat, p.value=p, #alternative="two-sided", 
    method="Difference-sign test of independence", data.name=dname,
    n=n, mu=mu, sigma=sigma)
  class(rval) <- "htest"
  rval
} #function


turningpoint.test <- function(x)
{
#check arguments
  if (!is.numeric(x)) #is x numeric?
    stop("x must be numeric.")
  if (NCOL(x)>1)
    stop("x is not a vector or univariate time series")
  dname <- deparse(substitute(x)) #name of data variable
  x <- c(x) #coerse to a vector and strip names

#Trim repeated values from the data
  x <- x[diff(x)!=0]

#Perform the test
  n <- length(x) #length of data set
  if (n<3)
    stop("there is not enough data to perform the turning point test")
  d1 <- diff(x[1:(n-1)]) #get first n-1 differences
  d2 <- diff(x[2:n]) #get last n-1 differences
  tp = sum((d1>0 & d2<0) | (d1<0 & d2>0)) #count turning points
  mu <- 2*(n-2)/3 #the expected number of turning points in n IID observations
  sigma <- sqrt((16*n-29)/90) #the variance of the number of turning points in n IID observations
  stat <- (tp-mu)/sigma #compute the test statistic which is asymptotically standard normal
  names(stat) <- "T" #name of statistic
  p <- 2*pnorm(abs(stat), lower.tail=FALSE) #get the p-value of the test statistic
#Pack return
  rval <- list(statistic=stat, p.value=p, #alternative="two-sided", 
    method="Turning point test of independence", data.name=dname,
    n=n, mu=mu, sigma=sigma)
  class(rval) <- "htest"
  rval
} #function


rank.test <- function(x)
{
#check arguments
  if (!is.numeric(x)) #is x numeric?
    stop("x must be numeric.")
  if (NCOL(x)>1)
    stop("x is not a vector or univariate time series")
  dname <- deparse(substitute(x)) #name of data variable
  x <- c(x) #coerse to a vector and strip names
#Trim repeated values from the data
  x <- x[diff(x)!=0]
#Perform the test
  n <- length(x) #length of data set
  if (n<2)
    stop("There is not enough data to perform a rank test.")
#  pairs <- 0
#  for (i in 1:(n-1))
#   pairs <- pairs + sum(x[i:(n-1)]<generate x[(i+1):n]) #count positive differences
#  pairs <- sum(sapply(2:n, function(i) sum(x[i:n]>x[i-1]))) #count all pairs of positive differences
pairs <- .C("CountIncreasingPairs", as.double(x), as.integer(n), double(n), count=double(1), overflow=integer(1))
if (pairs$overflow)
  if (pairs$count==2^32-1)
    warning("a 32-bit integer overflow has occurred and the results of rank.test should not be trusted")
  else
    warning("a 64-bit integer overflow has occurred and the results of rank.test should not be trusted")
  mu <- n*(n-1)/4 #the expected number of increasing pairs
  sigma <- sqrt(n*(n-1)*(2*n+5)/72) #the variance of the number of increasing pairs in n IID observations
  stat <- (pairs$count-mu)/sigma #compute the test statistic which is asymptotically standard normal
  names(stat) <- "R" #name of statistic
  p <- 2*pnorm(abs(stat), lower.tail=FALSE) #get the p-value of the test statistic
#Pack return
  rval <- list(statistic=stat, p.value=p, #alternative="two-sided", 
    method="Rank test of independence", data.name=dname,
    pairs=pairs$count, n=n, mu=mu, sigma=sigma)
  class(rval) <- "htest"
  rval
} #function


chisq.unif.test <- function(x, bins=NULL, interval=c(0,1), min.bin.size=10, all.inside=TRUE, rightmost.closed=TRUE, ...)
{
#check arguments
  if (!is.numeric(interval) ||  length(interval)!=2 || !all(is.finite(interval)) || any(is.na(interval)) || interval[1]>=interval[2])
    stop("interval must be a vector of 2 finite values whose first element is strictly smaller than its second element.")
  if (!is.numeric(x)) #is x numeric?
    stop("x must be numeric.")
  dname <- deparse(substitute(x)) #name of data variable
  if (!is.vector(x)) #is x a vector?
    x <- c(x) #no, coerse it to a vector
#  if (any(x<interval[1] | x>interval[2]))
#    stop(paste0("x must be a vector of values in the range ", interval[1], "-", interval[2], "."))
#Perform the test
  n <- length(x) #length of data set
  if (n<30)
    stop("there is not enough data.  At least 30 data points are needed")
#Declare function used for discretising the interval range
	discretise <- function(x, bins) findInterval(x, seq(interval[1], interval[2], length.out=bins+1), rightmost.closed=rightmost.closed, all.inside=all.inside)
#Discretise the specified interval 
  if (is.null(bins))
  {
    if (n>=200) bins <- 20
    else bins <- n/10
  	xdisc <- discretise(x, bins, ...)
  	xdisc <- xdisc[xdisc>=1 & xdisc<=bins]
  	counts <- table(xdisc)
  	while (bins>1 && min(counts)<min.bin.size)
  	{
    	bins <- bins-1
  		xdisc <- discretise(x, bins, ...)
  		xdisc <- xdisc[xdisc>=1 & xdisc<=bins]
  		counts <- table(xdisc)
  	} #while
  	if (bins==1)
  		stop("There is not enough data in each bin to perform the test.")
  }
  else
  {
  	xdisc <- discretise(x, bins, ...)
  	xdisc <- xdisc[xdisc>=1 & xdisc<=bins]
  	counts <- table(xdisc)
  } #if
#Perform test and pack return
  rval <- chisq.test(counts, ...)
  rval$method <- paste0("Discrete uniform(", format(interval[1], trim=TRUE), ",", format(interval[2], trim=TRUE), ") chi-squared test") 
  rval$parameter <- c(rval$parameter, a=interval[1], b=interval[2])
  rval$data.name <- dname
  rval$bins <- bins
  rval$min.bin.size <- min(counts)
  rval$interval <- interval
  rval
} #function

ks.unif.test <- function(x)
{
	rval <- ks.test(x, punif)
	rval$method <- "Uniform(0,1) Kolmogorov-Smirnov test"
	rval
} #function

lb.test <- function(x, ...) Box.test(x, type="Ljung-Box", ...)
