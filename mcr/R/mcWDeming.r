###############################################################################
##
## mcWDeming.r
##
## Function for computing weighted deming regression for two methods with  proportional errors.
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

#' Calculate Weighted Deming Regression
#'
#' Calculate weighted deming regression with iterative algorithm suggested by Linnet. 
#' This algorithm is avalaible only for positive values. But even in this case there is no guarantee that
#' the algorithm always converges. 
#'
#' @param X measurement values of reference method.
#' @param Y measurement values of test method.
#' @param error.ratio ratio between squared measurement errors of reference- and test method, necessary for Deming regression (Default is 1).
#' @param iter.max maximal number of iterations.
#' @param threshold threshold value.
#' @return a list with elements
#'  \item{b0}{intercept.}
#'  \item{b1}{slope.}
#'  \item{xw}{average of reference method values.}
#'  \item{iter}{number of iterations.}
#' @references  Linnet K.
#'              Evaluation of Regression Procedures for Methods Comparison Studies.
#'              CLIN. CHEM. 39/3, 424-432 (1993).
#' 
#'              Linnet K.
#'              Estimation of the Linear Relationship between the Measurements of two Methods with Proportional Errors.
#'              STATISTICS IN MEDICINE, Vol. 9, 1463-1473 (1990).
mc.wdemingConstCV <- function(X, Y, error.ratio, iter.max=30, threshold=0.000001)
{
  # Check validity of parameters
  
    stopifnot(is.numeric(X))
    stopifnot(is.numeric(Y))
    stopifnot(length(X) == length(Y))
    stopifnot(is.numeric(error.ratio))
    stopifnot(error.ratio > 0)
    stopifnot(is.numeric(iter.max))
    stopifnot(round(iter.max) == iter.max)
    stopifnot(iter.max > 0)
    stopifnot(is.numeric(threshold))
    stopifnot(threshold >= 0)
  
    # This algorithm often doesn't converge if there are negative
    # measurements in data set
  
    if (min(X)<0 | min(Y)<0)
    {
        return(paste("There are negative values in data set."))
    }
    else
    {
        # 1. Calculate  initial values
        #    (point estimates of unweighted deming regression)

        n <- length(X)
      	
      	mX <- mean(X)
      	mY <- mean(Y)
      	u <- sum((X-mX)^2)
      	q <- sum((Y-mY)^2)
      	p <- sum((X-mX)*(Y-mY))
      	
      	## initial values
      	
      	b1 <- ((error.ratio*q-u)+sqrt((u-error.ratio*q)^2+4*error.ratio*p^2))/(2*error.ratio*p)
      	b0 <- mean(Y)-b1*mean(X)
      	
        ## Iterative Algorithm
    
        i <- 0     # Number of iterations
        warn<-"no warnings"   # Warnings
      
        repeat
        {
            if (i >= iter.max) 
            {
                warning(paste("no konvergence after",iter.max,"iterations"))
                break
            }
        
            i<-i+1
          
            # Calculation of weights
            d <- Y-(b0+b1*X)
            XHAT <- X+(error.ratio*b1*d/(1+error.ratio*b1^2))
            YHAT <- Y-(d/(1+error.ratio*b1^2))
            W <- ((XHAT+error.ratio*YHAT)/(1+error.ratio))^(-2)
    
            # Calculation of regression coefficients
            XW <- sum(W*X)/sum(W)
            YW <- sum(W*Y)/sum(W)
            U <- sum(W*((X-XW)^2))
            Q <- sum(W*((Y-YW)^2))
            P <- sum(W*(X-XW)*(Y-YW))
              
            # Point estimates
            B1 <- (error.ratio*Q-U+sqrt((U-error.ratio*Q)^2+4*error.ratio*P^2))/(2*error.ratio*P)
            B0 <- YW-B1*XW
          
            # Stop condition
            if(abs(b1-B1) < threshold & abs(b0-B0) < threshold) 
                break
    
            # new values
            b1<-B1
            b0<-B0
        } # end of iterative algorithm
		
        return(list(b1=B1,b0=B0,iter=i,xw=XW, weight=W))
    }
}

