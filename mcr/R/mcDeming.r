###############################################################################
##
## mcDeming.R
##
## Function for computing Deming regression based method comparison.
##
## Copyright (C) 2011 Roche Diagnostics GmbH
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

#' Calculate Unweighted Deming Regression and Estimate Standard Errors
#' 
#' @param X measurement values of reference method.
#' @param Y measurement values of test method.
#' @param error.ratio ratio of measurement error of reference method to measurement error of test method.
#' @return a list with elements
#'  \item{b0}{intercept.}
#'  \item{b1}{slope.}
#'  \item{se.b0}{respective standard error of intercept.}
#'  \item{se.b1}{respective standard error of slope.}
#'  \item{xw}{average of reference method values.}
#' @references  Linnet K.
#'              Evaluation of Regression Procedures for Methods Comparison Studies.
#'              CLIN. CHEM. 39/3, 424-432 (1993).
#' 
#'              Linnet K.
#'              Estimation of the Linear Relationship between the Measurements of two Methods with Proportional Errors.
#'              STATISTICS IN MEDICINE, Vol. 9, 1463-1473 (1990).
mc.deming<-function(X, Y, error.ratio)
{    
    # Check validity of parameters
      
    stopifnot(!is.na(X))
    stopifnot(!is.na(Y))
    stopifnot(is.numeric(X))
    stopifnot(is.numeric(Y))
    stopifnot(length(X)==length(Y))
    stopifnot(length(X) > 0)
    stopifnot(!is.na(error.ratio))
    stopifnot(is.numeric(error.ratio))
    stopifnot(error.ratio>0)
    stopifnot(length(error.ratio) > 0)
    
    #--
    n <- length(X)
	mX <- mean(X)
	mY <- mean(Y)
	u <- sum((X-mX)^2)
	q <- sum((Y-mY)^2)
	p <- sum((X-mX)*(Y-mY))
	r <- p/sqrt(u*q)

	## Estimated points
	# [ Ref. K.Linnet. Estimation of the linear relationship between
    #        the measurements of two methods with  Proportional errors.
    #        STATISTICS IN MEDICINE, VOL. 9, 1463-1473 (1990)].
	#
	
	b1 <- ((error.ratio*q-u)+sqrt((u-error.ratio*q)^2+4*error.ratio*p^2))/(2*error.ratio*p)
	b0 <- mean(Y)-b1*mean(X)

    ## Standard error
    # [Ref. Strike, P. W. (1991) Statistical Methods in Laboratory Medicine.
    #       Butterworth-Heinemann, Oxford ].
  
    se.b1 <- sqrt(b1^2*(calcDiff(1,r^2)/r^2)/(n-2))
	se.b0 <- sqrt(se.b1^2*mean(X^2))

	return(list(b0=b0, b1=b1, se.b0=se.b0, se.b1=se.b1, xw=mX,  weight=rep(1,length(X))))
}
