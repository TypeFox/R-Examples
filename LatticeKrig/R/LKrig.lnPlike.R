# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2

LKrig.lnPlike <- function(Mc, Q, wy, residuals, weights, LKinfo) {
	# sanity check on lambda being ratio of sigma^2 and rho
	# if sigma and rho are passed.
	sigma <- LKinfo$sigma
	rho   <- LKinfo$rho
	lambda<- LKinfo$lambda
if (!is.na(rho)) {
		if ((sigma^2/rho) != lambda) {
			stop(" sigma, rho and lambda do not match up")
		}
	}
	#
	n <- nrow(wy)
	m <- dim(Q)[1]
	# find log determinant of reg matrix for use in the log likeihood
	lnDetReg <- 2 * sum(log(diag(Mc)))
	# log determinant of precision matrix.
	lnDetQ <- 2 * sum(log(diag(chol(Q, memory = LKinfo$choleskyMemory))))
	# now apply a miraculous determinant identity (Sylvester''s theorem)
	#  det( I + UV) = det( I + VU)    providing UV is square
# or more generally
#  det( A + UV) = det(A) det( I + V A^{-1}U)
#  or as we use it
#  ln(det( I + V A^{-1}U)) = ln( det(  A + UV)) - ln( det(A))
#
lnDetCov <- lnDetReg - lnDetQ + (n - m) * log(lambda) - sum(log(weights))
	# finding quadratic form
	# this uses a shortcut formula in terms of
    # the residuals	
###############DEBUG    
#c.mKrig <- weights * residuals/lambda
#y<- wy/ sqrt( weights)
#quad.form <- c(colSums(as.matrix(c.mKrig * y)))
#print( quad.form)
###############DEBUG
quad.form <- c(colSums(as.matrix(sqrt(weights)*wy*residuals/lambda) ) )
    # MLE estimate of rho and sigma
	# these are derived by assuming Y is  MN(  Ud, rho*M )
	rho.MLE <- quad.form/n
	shat.MLE <- sigma.MLE <- sqrt(lambda * rho.MLE)
	# the  log profile likehood with  rho.MLE  and  dhat substituted
	# leaving a profile for just lambda.
# note that this is _not_  -2*loglike just the log and
# includes the constants
lnProfileLike <- (-(n/2) - log(2 * pi) * (n/2) - (n/2) * log(rho.MLE) - 
		(1/2) * lnDetCov)
	# find log likelihood without profiling if sigma and rho have been passed.
	# NOTE: this assumes that  lambda == sigma^2/rho
if (!is.na(rho)) {
		lnLike <- (-(quad.form)/(2 * rho) - log(2 * pi) * (n/2) - (n/2) * 
			log(rho) - (1/2) * lnDetCov)
		lnLike.FULL <- sum(lnLike)
	} else {
		lnLike <- NA
		lnLike.FULL <- NA
	}
	# the full ln profile likelihood  is found 
	# by assuming the replicate fields (columns of y)
	# are independent and putting in a pooled MLE for rho.
rho.MLE.FULL <- mean(rho.MLE)
	sigma.MLE.FULL <- sqrt(lambda * rho.MLE.FULL)
	lnProfileLike.FULL <- sum(-(n/2) - log(2 * pi) * (n/2) - (n/2) * log(rho.MLE.FULL) - 
		(1/2) * lnDetCov)

	#
	return(list(lnProfileLike = lnProfileLike,
	                  rho.MLE = rho.MLE,
	                 shat.MLE = shat.MLE, 
	          		sigma.MLE = shat.MLE,
	          		    sigma = sigma,
	          		      rho = rho,
	          		   lnLike = lnLike,
	              lnLike.FULL = lnLike.FULL, 
           	     rho.MLE.FULL = rho.MLE.FULL,
           	   sigma.MLE.FULL = sigma.MLE.FULL,
           lnProfileLike.FULL = lnProfileLike.FULL, 
             		quad.form = quad.form,
             		 lnDetCov = lnDetCov))
}

