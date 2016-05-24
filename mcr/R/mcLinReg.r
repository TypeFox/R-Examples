###############################################################################
##
## mcLinReg.R
##
## Function for computing linear regression based method comparison.
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

#' Calculate ordinary linear Regression
#' and Estimate Standard Errors
#'
#' @param X measurement values of reference method.
#' @param Y measurement values of test method.
#' @return a list with elements
#'  \item{b0}{intercept.}
#'  \item{b1}{slope.}
#'  \item{se.b0}{respective standard error of intercept.}
#'  \item{se.b1}{respective standard error of slope.}
#'  \item{xw}{average of reference method values.}
#' @references  Neter J., Wassermann W., Kunter M.
#'              Applied Statistical Models. 
#'              Richard D. Irwing, INC., 1985.
mc.linreg <- function(X,Y) 
{
    ## Check validity of parameters
    stopifnot(!is.na(X))
    stopifnot(!is.na(Y))
    stopifnot(is.numeric(X))
    stopifnot(is.numeric(Y))
    stopifnot(length(X)==length(Y))
    stopifnot(length(X) > 0)

	## Number of data points
	n <- length(X)
    W <- rep(1,n)

    #--
  	XW <- sum(W*X)/sum(W)
    YW <- sum(W*Y)/sum(W)
    SXXW <- sum(W*X^2)-(sum(W*X)^2)/sum(W)
    SXYW <- sum(W*X*Y)-(sum(W*X)*sum(W*Y)/sum(W))

	## Point estimates

  	b1 <- SXYW/SXXW
  	b0 <- YW-b1*XW

	## Standard errors
	
    MSEW <- 1/(n-2)*sum(W*(Y-(b0+b1*X))^2)
    se.b1 <- sqrt(MSEW/SXXW)
  	se.b0 <- se.b1*sqrt(sum(W*X^2)/sum(W))

	## Return estimates and sd

    return(list(b0=b0,b1=b1,se.b0=se.b0,se.b1=se.b1,xw=XW, weight=rep(1,length(X))))
}

#' Calculate Weighted Ordinary Linear Regression
#' and Estimate Standard Errors
#'
#'  The weights of regression are taken as reverse squared values of the reference method, 
#'  that's why it is impossible to achieve the calculations for zero values. 
#'  
#' @param X measurement values of reference method.
#' @param Y measurement values of test method.
#' @return a list with elements.
#'  \item{b0}{intercept.}
#'  \item{b1}{slope.}
#'  \item{se.b0}{respective standard error of intercept.}
#'  \item{se.b1}{respective standard error of slope.}
#'  \item{xw}{weighted average of reference method values.}
#' @references  Neter J., Wassermann W., Kunter M.
#'              Applied Statistical Models. 
#'              Richard D. Irwing, INC., 1985.
mc.wlinreg <- function(X,Y) 
{
    ## Check validity of parameters
    stopifnot(!is.na(X))
    stopifnot(!is.na(Y))
    stopifnot(is.numeric(X))
    stopifnot(is.numeric(Y))
    stopifnot(length(X)==length(Y))
    stopifnot(length(X) > 0)

  	## Number of data points
	  n <- length(X)

    ## Weights
    W <- 1/X^2

    #--
  	XW <- sum(W*X)/sum(W)
    YW <- sum(W*Y)/sum(W)
    SXXW <- sum(W*X^2)-(sum(W*X)^2)/sum(W)
    SXYW <- sum(W*X*Y)-(sum(W*X)*sum(W*Y)/sum(W))

	## Point estimates

  	b1<-SXYW/SXXW
  	b0<-YW-b1*XW

	## Standard errors
	
    MSEW <- 1/(n-2)*sum(W*(Y-(b0+b1*X))^2)
    se.b1 <- sqrt(MSEW/SXXW)
  	se.b0 <- se.b1*sqrt(sum(W*X^2)/sum(W))

	## Return estimates and sd

    return(list(b0=b0,b1=b1,se.b0=se.b0,se.b1=se.b1,xw=XW, weight=W))
}

