#******************************************************************************* 
#
# Bayesian Regression and Adaptive Sampling with Gaussian Process Trees
# Copyright (C) 2005, University of California
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (rbgramacy@ams.ucsc.edu)
#
#*******************************************************************************


## default.itemps:
##
## create a default inverse temperature ladder for importance
## tempering (IT) together with pseudo-prior and parameters
## for calibrating it by stochastic approximation.  There are three
## choices of ladder as specified by type

"default.itemps" <-
function(m=40, type=c("geometric", "harmonic", "sigmoidal"),
         k.min=0.1, c0n0=c(100,1000), lambda=c("opt", "naive", "st"))
{
  ## check m argument
  if(length(m) != 1 || m <= 0)
    stop("m should be a positive integer")
  
  ## check type argument
  type <- match.arg(type)

  ## check k.min argument
  if(length(k.min) != 1 || k.min >= 1 || k.min < 0)
    stop("k.min should be a integer satisfying 0 <= k.min < 1")

  ## check the c0n0 argument
  if(length(c0n0) != 2 || !prod(c0n0 >= 0))
    stop("c0n0 should be a nonnegative 2-vector")

  ## check the lambda argument
  lambda <- match.arg(lambda)

  ## check if importance sampling only
  if(m == 1) return(list(c0n0=c(0,0), k=k.min, pk=1, lambda="naive"))
  
  if(type == "geometric") {

    ## calculate the delta for the geometric which reaches
    ## k.min in m steps
    delta <- k.min^(1/(1-m)) - 1

    ## geometric temperature ladder
    i <- 1:m
    k <- (1+delta)^(1-i)
    
  } else if(type == "harmonic") {

    ## calculate the delta for the geometric which reaches
    ## k.min in m steps
    delta <- ((1/k.min) - 1)/(m-1)

    ## harmonic temperature ladder
    i <- 1:m
    k <- 1/(1+ delta*(i-1))

  } else { ## sigmoid

    ## calculate the indices which provide the sigmoid which
    ## begins at 1 and ends at k.min with m steps
    x <- c(1,k.min)
    ends <- log((1.01-x)/x)
    t <- seq(ends[1], ends[2], length=m)

    ## logistic/sigmoid temperature ladder
    k <- 1.01 - 1.01/(1+exp(-t))
    
  }

  ## return the generated ladder, as above, with a vector of
  ## observation counts for tgp to update
  return(list(c0n0=c0n0, k=k, pk=rep(1/m, m), lambda=lambda))
}


## check.itemps:
##
## check the itemps create by hand or from default.itemps or
## as modified by tgp or predict tgp and assembled inside
## the tgp.postprocess function

"check.itemps" <- 
function(itemps, params)
{
  ## if null, then just make one temperature (1.0) with all the prob
  if(is.null(itemps)) return(c(1,0,0,1,1,0,1))

  ## if it is a list or a data frame
  else if(is.list(itemps) || is.data.frame(itemps)) {
    
    ## get the four fields
    c0n0 <- itemps$c0n0
    pk <- itemps$pk
    lambda <- itemps$lambda
    k <- itemps$k
    counts <- itemps$counts
    
    ## check for non-null k
    m <- length(k)
    if(m == 0) stop("must specify k vector in list")
    
    ## check for null pk
    if(is.null(pk)) pk <- rep(1/m, m)
    
    ## check the dims are right
    if(m != length(pk))
      stop("length(itemps$k) != length(itemps$pk)")
    
    ## put into decreasing order
    o <- order(k, decreasing=TRUE)
    k <- k[o]
    pk <- pk[o]
    
    ## checks k
    if(prod(k >= 0)!=1) stop("should have 0 <= itemps$k")
    if((m > 1 || k != 1) && params$bprior != "b0")
      warning("recommend params$bprior == \"b0\" for itemps$k != 1",
               immediate.=TRUE)
  
    ## checks for pk
    if(prod(pk > 0)!=1) stop("all itemps$pk should be positive")

    ## init and checks for c0n0
    if(! is.null(c0n0)) {
      if(length(c0n0) != 2 || !prod(c0n0 >= 0))
        stop("itemps$c0n0 should be a nonnegative 2-vector")
    } else c0n0 <- c(100,1000)

    ## check lambda
    if(! is.null(lambda)) {
      if(lambda == "opt") lambda <- 1
      else if(lambda == "naive") lambda <- 2
      else if(lambda == "st") {
        if(k[1] != 1.0) stop("cannot use lambda=\"st\" when itemps$k[1] != 1.0\n")
         lambda <- 3
      } else stop(paste("lambda = ", lambda, "is not valid\n", sep=""))
    } else lambda <- 1 

    ## check the counts vector
    if(! is.null(counts)) {
      if(m != length(counts)) stop("length(itemps$k) != length(itemps$counts)")
    } else counts <- rep(0,m)
    
    ## return a double-version of the ladder
    return(c(m, c0n0, k, pk, counts, lambda))
  }

  ## if it is a matrix
  else if(is.matrix(itemps)) {

    ## check dims of matrix
    if(ncol(itemps) != 2) stop("ncol(itemps) should be 2")

    ## get the two fields
    pk <- itemps[,2]
    k <- itemps[,1]
    m <- length(k)

    ## put into decreasing order
    o <- order(k, decreasing=TRUE)
    k <- k[o]
    pk <- pk[o]
    
    ## checks k
    if(prod(k >= 0)!=1) stop("should have 0 <= itemps[,1]")
    if((m > 1 || k != 1) && params$bprior != "b0")
      warning("recommend params$bprior == \"b0\" for itemps[,1] != 1",
              immediate.=TRUE)

    ## checks for pk
    if(prod(pk > 0)!=1) stop("all probs in itemps[,2] should be positive")

    ## return a double-version with a counts vector at the end
    return(c(m, 100, 1000, k, pk, 1, rep(0,m)))
  }

  ## if itemps is a vector
  else if(is.vector(itemps)) {

    ## get length of inverse temperature ladder
    m <- length(itemps)
    
    ## checks for itemps
    if(prod(itemps >= 0)!=1) stop("should have 0 <= itemps ")
    if((length(itemps) > 1 || itemps != 1) && params$bprior != "b0")
      warning("recommend params$bprior == \"b0\" for itemps != 1",
               immediate.=TRUE)

    ## return a double-version with a counts vector at the end
    return(c(m, 100, 1000, itemps, rep(1/m, m), 1, rep(0,m)))
  }
  else stop("invalid form for itemps")
}


## hist2bar:
##
## make a barplot to compare the (discrete) histograms
## of each column of the input argument x

hist2bar <- function(x)
{
  ## make a matrix
  if(is.vector(x)) x <- matrix(x, ncol=1)
  
  ## calculate the number of, and allocate the space for,
  ## the bins, b, of the histogram
  r <- range(as.numeric(x))
  b <- matrix(0, ncol=ncol(x), nrow=r[2]-r[1]+1)

  ## calculate the histogram height of each bin
  for(i in r[1]:r[2])
    for(j in 1:ncol(x))
      b[i-r[1]+1,j] <- sum(x[,j] == i)

  ## make have thr right data.frame format so that
  ## it will place nice with the barplot function,
  ## and return
  b <- data.frame(b)
  row.names(b) <- r[1]:r[2]
  return(t(b))
}


## itemps.barplot:
##
## make a histogram (via barplot) of the number of times
## each inverse-temperature was visited in the ST-MCMC
## chain.  Requires that traces were collected

itemps.barplot <- function(obj, main=NULL, xlab="itemps",
                           ylab="counts", plot.it=TRUE, ...)
{
  ## check to make sure traces were collected
  if(is.null(obj$trace)) 
    stop(paste("no traces in tgp-object;", 
               "re-run the b* function with argument \"trace=TRUE\""))

  ## check to make sure tempering was used
  if(is.null(obj$itemps)) stop("no itemps in tgp-object")

  ## create a bin for each inverse-temperature
  bins <- rep(0,length(obj$itemps$k))

  ## count and store the number in the first bin
  m <- obj$trace$post$itemp == obj$itemps$k[1]
  bins[1] <- sum(m)

  ## count and store the number in the rest of the bins
  for(i in 2:length(obj$itemps$k)) {
    m <- obj$trace$post$itemp == obj$itemps$k[i]
    if(sum(m) == 0) next;
    bins[i] <- sum(m)
  }

  ## make into a data frame for convenient barplotting
  bins <- data.frame(bins)
  row.names(bins) <- signif(obj$itemps$k,3)
  
  ## make the barplot histogram
  if(plot.it==TRUE) {
    smain <- paste(main, "itemp counts")
    barplot(t(bins), xlab=xlab, ylab=ylab, ...)
  }

  ## return the barplot structure for plotting later
  return(invisible(bins))
}
