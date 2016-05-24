#*******************************************************************************
#
# Local Approximate Gaussian Process Regression
# Copyright (C) 2013, The University of Chicago
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
# Questions? Contact Robert B. Gramacy (rbgramacy@chicagobooth.edu)
#
#*******************************************************************************


## approxGP.seq:
##
## fit a sequence of approx GPs, following the method progression specified
## -- assumes 1/1000 nugget

aGP.seq <- function(X, Z, XX, d, methods=rep("alc", 2), M=NULL, ncalib=0, 
	## aGP=aGP, 
	...)
  {
  	## base case for not emulating, M function must be defined
  	if(length(methods) == 1 && methods == FALSE) {
  		nc <- ncol(X)
  		return(data.frame(mean=M(XX[,1:(nc-ncalib)], XX[,(nc-ncalib+1):nc]), var=0))
  	}

  	M.hat <- NULL
	for(m in methods) {
		M.hat <- aGP(X=X, Z=Z, XX=XX, method=m, d=d, ...) 
		if(length(d$start) == 1) d$start <- rep(d$start, nrow(XX))
		nna <- is.finite(M.hat$mle$d)
		d$start[nna] <- M.hat$mle$d[nna]
	}
	return(M.hat)
  }


## discrep.est:
##
## switch function for estimating a GP discrepency, or simply
## using zero-mean IID and estimating variance (bias=FALSE)

discrep.est <- function(X, Y, Yhat, d, g, bias=TRUE, clean=TRUE)
	{
		## WHEN REPS > 1 COULD USE DUPLICATES TO GET INITIAL NUGGET ESTIMATE

		## remove columns from X which have fewer than three
		lun <- apply(X, 2, function(x) { length(unique(x)) })
		lun2 <- which(lun <= 2)
		if(length(lun2) > 0) X <- X[,-lun2]

		## calculate discrepency: possibly dealing with replicates
		reps <- length(Y)/nrow(X)
		if(reps != 1) {
			X <- matrix(rep(t(X), reps), byrow=TRUE, ncol=ncol(X))
			Ym <- Y - rep(Yhat, reps)
		} else 
			Ym <- Y - Yhat

		## estimate discrepency or not
		if(bias) { 
			gp <- newGP(X, Ym, d=d$start, g=g$start, dK=TRUE)
			cmle <- jmleGP(gp, drange=c(d$min, d$max), grange=c(g$min, g$max), 
			 			   dab=d$ab, gab=g$ab)
		} else { ## case for no discrepency estimated - just variance
			gp <- newGP(X, Ym, d=0, g=0)
			cmle <- list(d=0, g=0, dits=0, gits=0)		
		}

		## calculate the log likelihood
		cmle$ll <- llikGP(gp, dab=d$ab, gab=g$ab)
		
		## clean up
		if(clean) deleteGP(gp)
		else cmle$gp <- gp

		## done
		return(cmle)
	}


## fcalib:
##
## function evaluation by an optimizer in order to find u.hat
## for calibration -- returns the negative log likelihood

fcalib <- function(u, XU, Z, X, Y, da, d, g, uprior=NULL,
				   methods=rep("alc", 2), M=NULL, bias=TRUE,
				   omp.threads=1, save.global=FALSE, verb=1)
  {
	## fit an approxGP emulator to the computer model data
	Xu <- cbind(as.matrix(X), matrix(rep(as.numeric(u), nrow(X)), 
		ncol=length(u), byrow=TRUE))

	## loop over methods
	M.hat <- aGP.seq(X=XU, Z=Z, XX=Xu, d=da, methods=methods, M=M, 
		ncalib=ncol(XU)-ncol(X), omp.threads=omp.threads, verb=0)

	## find mle parameters for discrepency model
	cmle <- discrep.est(X, Y, M.hat$mean, d, g, bias, TRUE)

	## possibly add in prior for u
	if(!is.null(uprior)) cmle$ll <- cmle$ll + uprior(u)
	
	## progress meter
	if(verb > 0) 
		cat("u=(", paste(u, collapse=", "), "), its=", cmle$dits+cmle$gits, 
			" (d,g)=(", cmle$d, ",", cmle$g, "), llik=", cmle$ll, "\n", sep="")

	## possibly save to a global variable
	if(!is.logical(save.global)) {
		if(exists("fcalib.save", envir=save.global)) {
			assign("fcalib.save", rbind(get("fcalib.save", save.global), 
				c(u, cmle$ll)), save.global)
		} else {
			assign("fcalib.save", data.frame(u=matrix(u, nrow=1), ll=cmle$ll), 
				save.global)
		}
	} else if(is.logical(save.global) && save.global == TRUE)
		stop("when not FALSE, save.global must be an environment")

	## done
	return(-cmle$ll)
  }
