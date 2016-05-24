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


## tgp.design:
##
## choose howmany of Xcand candidate locations according to a treed
## D-optimal design using the MAP tree contained in the tgp-class
## object.  iter specifies the number of iterations in the stochastic
## ascent method

"tgp.design" <-
function(howmany, Xcand, out, iter=5000, verb=0)
{
  ## get partitioned candidates and dat locaitons
  Xcand.parts <- partition(Xcand, out)
  X.parts <- partition(out$X, out)

  ## initialize selected candidates to none
  XX <- NULL
  
  ## subsample some from each partition
  cat(paste("\nsequential treed D-Optimal design in ", 
            length(Xcand.parts), " partitions\n", sep=""))
  for(i in 1:length(Xcand.parts)) {
    nn <- ceiling(howmany*(nrow(Xcand.parts[[i]]))/(nrow(Xcand)))
    if(verb > 0)
      cat(paste("dopt.gp (", i, ") choosing ", nn, " new inputs from ", 
                nrow(Xcand.parts[[i]]), " candidates\n", sep=""))
    dout <- dopt.gp(nn, X.parts[[i]], Xcand.parts[[i]], iter, max(verb-1,0));
    XX <- rbind(XX, dout$XX)
  }
  
  return(XX)
}


## tgp.partition:
##
## group X into a list containg each region as partitioned
## by the tree -- i is used to index root of the tree structure
## when applied recursively.  This functio is used exclusively
## by the partition function below

"tgp.partition" <-
function(X, tree, i)
{
  ## error or leaf node
  if(length(X) == 0) { stop("no X's found in partition\n") }
  if(tree$var[i] == "<leaf>") return(list(X));
  
  ## make sure X is a matrix
  if(is.null(nrow(X))) X <- matrix(X, ncol=1)
  
  ## gather the appropriate operations from the ith tree node
  var <- as.integer(as.character(tree$var[i]))+1
  gt <- (1:nrow(X))[X[,var] > tree$val[i]]
  leq <- setdiff(1:nrow(X), gt)
  
  ## calculate the left and right tree node rows
  l <- (1:nrow(tree))[tree$rows == 2*tree$rows[i]]
  r <- (1:nrow(tree))[tree$rows == 2*tree$rows[i]+1]
  
  ## recurse on left and right subtrees
  if(length(leq) > 0) Xl <- tgp.partition(as.matrix(X[leq,]), tree, l)
  else Xl <- NULL
  if(length(gt) > 0) Xr <- tgp.partition(as.matrix(X[gt,]), tree, r)
  else Xr <- NULL

  return(c(Xl,Xr))
}


## partition:
##
## return a list of X location in each region of the MAP
## treed partition contained in the tgp-class object in out

"partition" <-
function(X, out)
{
  m <- which.max(out$posts$lpost)
  tree <- out$trees[[out$posts$height[m]]]
  
  return(tgp.partition(X, tree, 1))
}


## dopt.gp:
##
## create a sequential D-optimal design of size under a GP model
## from candidates Xcand assuming that X locations are already in
## the design.  The stochastic ascent algorithm uses iter rounds.
## Uses a C-side routine via .C

"dopt.gp" <-
function(nn, X=NULL, Xcand, iter=5000, verb=0)
{
  if(nn == 0) return(NULL);

  ## check iterations
  if(length(iter) != 1 && iter <= 0)
    stop("iter must be a positive integer")

  ## check Kverbiterations
  if(length(verb) != 1 && iter < 0)
    stop("verb must be a non-negative integer")
  
  ## check X inputs
  Xnames <- names(X)
  X <- check.matrix(X)$X

  ## check the Xcand inputs
  if(is.null(Xcand)) stop("XX cannot be NULL")
  Xcand <- check.matrix(Xcand)$X

  ## check if X is NULL 
  if(!is.null(X)) {
    n <- nrow(X); m <- ncol(X)
    X <- t(X) ## for row-major in .C
  } else { n <- 0; m <- ncol(Xcand) }

  ## check that cols of Xcand match X
  if(ncol(Xcand) != m) stop("mismatched column dimension of X and Xcand");
  ncand <- nrow(Xcand)

  ## reduce nn if it is too big
  if(nn > nrow(Xcand)) {
    warning("nn greater than dim(Xcand)[1]");
    nn <- nrow(Xcand);
  }

  ## choose a random state for the C code
  state <- sample(seq(0,999), 3)

  ## run the C code
  ll <- .C("dopt_gp", 
           state = as.integer(state),
           nn = as.integer(nn),
           ## transpose of X is taken above
           X = as.double(X),
           n = as.integer(n),
           m = as.integer(m),
           Xcand = as.double(t(Xcand)),
           ncand = as.integer(ncand),
           iter = as.integer(iter),
           verb = as.integer(verb),
           fi = integer(nn),
           PACKAGE="tgp"
           )
  
  ## deal with X, and names of X
  ll$X <- framify.X(ll$X, Xnames, m)
  ll$Xcand <- framify.X(ll$Xcand, Xnames, m)
  ll$XX <- ll$Xcand[ll$fi,]
  if(is.matrix(Xcand)) ll$XX <- matrix(as.matrix(ll$XX), ncol=ncol(Xcand))

  ## dont return some of the things used by C
  ll$n <- NULL; ll$m <- NULL; ll$state <- NULL
  
  return(ll)
}

