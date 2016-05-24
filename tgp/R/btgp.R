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


## btgp:
##
## tgp implementation of the Bayesian treed Gaussian process model

"btgp" <-
function(X, Z, XX=NULL,
        meanfn="linear", bprior="bflat", corr="expsep", tree=c(0.5,2), 
	BTE=c(2000,7000,2), R=1, m0r1=TRUE, linburn=FALSE, itemps=NULL,
	pred.n=TRUE, krige=TRUE, zcov=FALSE, Ds2x=FALSE, improv=FALSE, 
        sens.p=NULL, nu=1.5, trace=FALSE, verb=1, ...)
{
  n <- nrow(X)
  if(is.null(n)) { n <- length(X); X <- matrix(X, nrow=n); d <- 1 }
  else { d <- ncol(X) }
  params <- tgp.default.params(d, meanfn=meanfn, corr=corr, ...)
  params$bprior <- bprior
  params$tree[1:length(tree)] <- tree
  params$gamma <- c(0,0.2,0.7)	# no llm
  if(corr == "matern") params$nu<-nu
  if(linburn && corr == "sim") stop("cannot do linburn for SIM model")
  return(tgp(X,Z,XX,BTE,R,m0r1,linburn,params,itemps,pred.n,krige,zcov, 
             Ds2x,improv,sens.p,trace,verb))
}


## bcart:
##
## tgp implementation of the Bayesian CART model of Chipman et al

"bcart" <-
function(X, Z, XX=NULL, bprior="bflat", tree=c(0.5,2), BTE=c(2000,7000,2), 
	R=1, m0r1=TRUE, itemps=NULL, pred.n=TRUE, krige=TRUE, zcov=FALSE, 
	Ds2x=FALSE, improv=FALSE, sens.p=NULL, trace=FALSE, verb=1, ...)
{
  return(btlm(X,Z,XX,meanfn="constant", bprior,tree,BTE,R,m0r1,itemps,pred.n,krige,
              zcov,Ds2x,improv,sens.p,trace,verb,...))
}


## bgp:
##
## tgp implementation of a Bayesian Gaussian process model

"bgp" <-
function(X, Z, XX=NULL, meanfn="linear", bprior="bflat", corr="expsep",
         BTE=c(1000,4000,2), R=1, m0r1=TRUE, itemps=NULL, pred.n=TRUE, 
         krige=TRUE, zcov=FALSE, Ds2x=FALSE, improv=FALSE, sens.p=NULL, nu=1.5,
         trace=FALSE, verb=1,  ... )
{
  n <- dim(X)[1]
  if(is.null(n)) { n <- length(X); X <- matrix(X, nrow=n); d <- 1 }
  else { d <- dim(X)[2] }
  params <- tgp.default.params(d, meanfn=meanfn, corr=corr,...)
  params$bprior <- bprior
  params$tree[1:2] <- c(0,0) # no tree
  params$gamma <- c(0,0.2,0.7)	# no llm
  if(corr == "matern") params$nu <- nu
  return(tgp(X,Z,XX,BTE,R,m0r1,FALSE,params,itemps,pred.n,krige,zcov,Ds2x,
             improv,sens.p,trace,verb))
}


## bgpllm:
##
## tgp implementation of a Bayesian Gaussian Process with
## jumps to the Limiting Linear Model

"bgpllm" <-
function(X, Z, XX=NULL, meanfn="linear", bprior="bflat", corr="expsep",
         gamma=c(10,0.2,0.7), BTE=c(1000,4000,2), R=1, m0r1=TRUE,
         itemps=NULL, pred.n=TRUE, krige=TRUE, zcov=FALSE,
	 Ds2x=FALSE, improv=FALSE, sens.p=NULL, nu=1.5, trace=FALSE, verb=1, ...)
{
  n <- dim(X)[1]
  if(is.null(n)) { n <- length(X); X <- matrix(X, nrow=n); d <- 1 }
  else { d <- dim(X)[2] }
  params <- tgp.default.params(d, meanfn=meanfn, corr=corr, ...)
  params$bprior <- bprior
  params$gamma <- gamma
  params$tree[1:2] <- c(0,0) # no tree
  if(corr == "matern"){ params$nu <- nu; }
  if(corr == "mrexpsep"){ stop("Sorry, the LLM is not yet available for corr=\"mrexpsep\"")}
  if(corr == "sim"){ stop("Sorry, the LLM is not available for corr=\"sim\"")}
  return(tgp(X,Z,XX,BTE,R,m0r1,FALSE,params,itemps,pred.n,krige,zcov,Ds2x,
             improv,sens.p,trace, verb))
}


## blm:
##
## tgp implementation of a Bayesian hierarchical Linear Model

"blm" <-
function(X, Z, XX=NULL, meanfn="linear", bprior="bflat",
         BTE=c(1000,4000,3), R=1, m0r1=TRUE, itemps=NULL, pred.n=TRUE, 
         krige=TRUE, zcov=FALSE, Ds2x=FALSE, improv=FALSE, sens.p=NULL, trace=FALSE, 
	 verb=1, ...)
{
  n <- dim(X)[1]
  if(is.null(n)) { n <- length(X); X <- matrix(X, nrow=n); d <- 1 }
  else { d <- dim(X)[2] }
  params <- tgp.default.params(d, meanfn=meanfn, ...)
  params$bprior <- bprior
  params$tree[1:2] <- c(0,0) # no tree
  params$gamma <- c(-1,0.2,0.7)	# force llm
  params$nug.p <- 0 ## force a fixed nugget
  params$gd[1] <- 0 ## fix the nugget at zero
  return(tgp(X,Z,XX,BTE,R,m0r1,FALSE,params,itemps,pred.n,
             krige,zcov,Ds2x,improv,sens.p,trace,verb))
}


## btgpllm:
##
## tgp implementation of a Bayesian treed Gaussian Process model
## with jumps to the Limiting Linear Model

"btgpllm" <-
function(X, Z, XX=NULL, meanfn="linear", bprior="bflat", corr="expsep",
         tree=c(0.5,2), gamma=c(10,0.2,0.7), BTE=c(2000,7000,2), R=1, m0r1=TRUE,
         linburn=FALSE, itemps=NULL, pred.n=TRUE, krige=TRUE, zcov=FALSE, Ds2x=FALSE, 
	 improv=FALSE, sens.p=NULL, nu=1.5, trace=FALSE, verb=1, ...)
{
  n <- nrow(X)
  if(is.null(n)) { n <- length(X); X <- matrix(X, nrow=n); d <- 1 }
  else { d <- ncol(X) }
  params <- tgp.default.params(d, meanfn=meanfn, corr=corr,...)
  params$bprior <- bprior
  params$tree[1:length(tree)] <- tree
  params$gamma <- gamma
  if(corr == "matern"){ params$nu <- nu }
  if(corr == "mrexpsep"){ stop("Sorry, the LLM  is not yet available for corr=\"mrexpsep\"")}
  if(corr == "sim"){ stop("Sorry, the LLM is not available for corr=\"sim\"")}
  return(tgp(X,Z,XX,BTE,R,m0r1,linburn,params,itemps,pred.n,krige,zcov,
             Ds2x,improv,sens.p,trace,verb))
}


"btlm" <-
function(X, Z, XX=NULL, meanfn="linear", bprior="bflat",
         tree=c(0.5,2), BTE=c(2000,7000,2), R=1, m0r1=TRUE, itemps=NULL, 
         pred.n=TRUE, krige=TRUE, zcov=FALSE, Ds2x=FALSE, improv=FALSE, 
         sens.p=NULL, trace=FALSE, verb=1, ...)
{
  n <- nrow(X)
  if(is.null(n)) { n <- length(X); X <- matrix(X, nrow=n); d <- 1 }
  else { d <- ncol(X) }
  params <- tgp.default.params(d, meanfn=meanfn, ...)
  params$bprior <- bprior
  params$tree[1:length(tree)] <- tree
  params$gamma <- c(-1,0.2,0.7)	# no llm
  params$nug.p <- 0 ## force a nugget 
  params$gd[1] <- 0 ## fix the nugget at zero
  return(tgp(X,Z,XX,BTE,R,m0r1,FALSE,params,itemps,pred.n,krige,zcov,
             Ds2x,improv,sens.p, trace,verb))
}

