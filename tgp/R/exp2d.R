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

## for R CMD CHECK
if(getRversion() >= "2.15.1") utils::globalVariables("exp2d", package="tgp")


## exp2d.Z:
##
## sample from he 2-d exponential data at locations X with
## normal mean-zero random deviates with sd specified

"exp2d.Z" <-
function(X, sd=0.001)
{
if(is.null(X)) return(NULL);
if(is.null(ncol(X))) X <- matrix(X, ncol=length(X))

## check the number of columns
if(ncol(X) != 2)
stop(paste("X should be a matrix (or data frame) with 2 columns, you have",
       ncol(X)))

## calculate the Z data
Ztrue <- X[,1] * exp(- X[,1]^2 - X[,2]^2)

## add randomness for random sample
Z <- Ztrue + rnorm(nrow(X),mean=0,sd=sd)

## return a data frame object
return(data.frame(Z=Z,Ztrue=Ztrue))
}


## exp2d.rand:
##
## samplig from the 2-d exponential data using the data file
## or a d-optimal (and/or) LH design and random evaluations from
## the exp2d.Z function.  n1 samples are taken from the interesting
## region, and n2 from outside.

"exp2d.rand" <-
function(n1=50, n2=30, lh=NULL, dopt=1)
{
## check the sanity of the inputs
if(n1 < 0 || n2 < 0) { stop("n1 and n2 must be >= 0") }

## use Latin Hybpercube sampling
if(!is.null(lh)) {

## start with the interesting region
Xcand <- lhs(n1*dopt, rbind(c(-2,2), c(-2,2)))
if(dopt > 2) { X <- dopt.gp(n1, NULL, Xcand)$XX }
else { X <- Xcand }

## check if n2 is a 1-vector or a 3-vector
if(length(n2) == 1) n2 <- rep(ceiling(n2/3), 3)
else if(length(n2 != 3))
stop(paste("length of n2 should be 1 or 3, you have",
	 length(n2)))

## check validity of dopt
if(length(dopt) != 1 || dopt < 1)
stop(paste("dopt should be a scalar >= 1, you have", dopt))

## do the remaining three (uninteresting) quadtants
Xcand <- lhs(n2[1]*dopt, rbind(c(2,6), c(-2,2)))
Xcand <- rbind(Xcand, lhs(n2[2]*dopt, rbind(c(2,6), c(2,6))))
Xcand <- rbind(Xcand, lhs(n2[3]*dopt, rbind(c(-2,2), c(2,6))))

## see if we need d-optimal subsample
if(dopt > 2) { X <- rbind(X, dopt.gp(sum(n2), NULL, Xcand)$XX) }
else { X <- rbind(X, Xcand) }

## calculate the Z data
Zdata <- exp2d.Z(X);
Ztrue <- Zdata$Ztrue; Z <- Zdata$Z

## now get the size of the XX vector (for each quadtant)

if(length(lh) == 1) lh <- rep(ceiling(lh/4), 4)
else if(length(lh) != 4)
stop(paste("length of lh should be 0 (for grid), 1 or 4, you have",
	 length(lh)))

## fill the XX vector
XX <- lhs(lh[1]*dopt, rbind(c(-2,2), c(-2,2)))
XX <- rbind(XX, lhs(lh[2]*dopt, rbind(c(2,6), c(-2,2))))
XX <- rbind(XX, lhs(lh[3]*dopt, rbind(c(2,6), c(2,6))))
XX <- rbind(XX, lhs(lh[4]*dopt, rbind(c(-2,2), c(2,6))))

## see if we need d-optimal subsample
if(length(X) > 0 && dopt > 2) {
XX <- dopt.gp(sum(lh), X, XX)$XX
}

## calculate the ZZ data
ZZdata <- exp2d.Z(XX);
ZZtrue <- ZZdata$Ztrue; ZZ <- Zdata$Z

} else {

## make sure we have enough data to fulfill the request
if(n1 + n2 >= 441) { stop("n1 + n2 must be <= 441") }

## dopt = TRUE doesn't make sense here
if(dopt != 1) { warning("argument dopt != 1 only makes sens when !is.null(lh)") }

## load the data
data(exp2d, envir=environment()); n <- dim(exp2d)[1]

    ## get the X columns
    si <- (1:n)[1==apply(exp2d[,1:2] <= 2, 1, prod)]
    s <- c(sample(si, size=n1, replace=FALSE), 
           sample(setdiff(1:n, si), n2, replace=FALSE))
    X <- as.matrix(exp2d[s,1:2]);

    ## get the XX predictive columns
    ss <- setdiff(1:n, s)
    XX <- exp2d[ss, 1:2];

    ## read the Z response columns
    Z <- as.vector(exp2d[s,3]);
    Ztrue <- as.vector(exp2d[s,4]);

    ## read the ZZ response columns
    ZZ <- as.vector(exp2d[ss,3]);
    ZZtrue <- as.vector(exp2d[ss,4]);
  }

  return(list(X=X, Z=Z, Ztrue=Ztrue, XX=XX, ZZ=ZZ, ZZtrue=ZZtrue))
}
