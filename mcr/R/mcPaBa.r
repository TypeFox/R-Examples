###############################################################################
##
## mcPaBa.R
##
## Functions for computing Passing Bablok regression.
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

#' Calculate Matrix of All Pair-wise Slope Angles
#' 
#' This is a very slow R version. It should not be called except for debugging purposes. 
#'
#' @param X measurement values of reference method.
#' @param Y measurement values of test method.
#' @param posCor should the algorithm assume positive correlation, i.e. symmetry around slope 1? 
#' @return Upper triangular matrix of slopes for all point combinations. Slopes in radian.
mc.calcAngleMat.R <- function(X, Y, posCor=TRUE)
{
	## Check validity of parameters
	stopifnot(is.numeric(X))
    stopifnot(is.numeric(Y))
    stopifnot(length(X)==length(Y))
	
	## Set up
	nData <- length(X)
	angM <- matrix(NA,nrow=nData,ncol=nData)
	
	## Calculate angles
	for(j in 1:(nData-1)) 
    {
		for(k in (j+1):nData) 
        {
			## Calculate x and y difference 
			dx <- calcDiff(X[k],X[j])
			dy <- calcDiff(Y[k],Y[j])	
			if(dx!=0) {
			## x and y != 0 
				angM[j,k] <- atan(dy/dx)
			}
			else if(dy!=0) 
            {
			## x==0, y!=0
			## only positive infinity for pos correlated
				if(posCor) 
                    angM[j,k] <- pi/2
				else 
                    angM[j,k] <- -pi/2
			}
			## else dx == dy == 0 => leave angM as NA, slope undefined 
		}
	}
	return(angM)
}

#' Calculate Matrix of All Pair-wise Slope Angles
#' 
#' This version is implemented in C for computational efficiency.
#'
#' @param X measurement values of reference method.
#' @param Y measurement values of test method.
#' @param posCor should algorithm assume positive correlation, i.e. symmetry around slope 1? 
#' @return Upper triangular matrix of slopes for all point combinations. Slopes in radian.
mc.calcAngleMat <- function(X,Y,posCor=TRUE) 
{
    ## Check validity of parameters
    stopifnot(is.numeric(X))
    stopifnot(is.numeric(Y))
    stopifnot(length(X)==length(Y))
    stopifnot(!is.na(posCor))
    ## Call C function
	ans <- .Call("calcAngleMat",X,Y,posCor)
	return(ans)
}


#' Passing-Bablok Regression
#'
#' @param angM upper triangular matrix of slopes for all point combinations. Slopes in radian.
#' @param X measurement values of reference method
#' @param Y measurement values of test method
#' @param alpha numeric value specifying the 100(1-alpha)% confidence level
#' @param posCor should algorithm assume positive correlation, i.e. symmetry around slope 1?
#' @param calcCI should confidence intervals be computed?
#' @param slope.measure angular measure of pairwise slopes  (see \code{\link{mcreg}} for details).\cr   
#'          \code{"radian"} - for data sets with even sample numbers median slope is calculated as average of two central slope angles.\cr
#'          \code{"tangent"} - for data sets with even sample numbers median slope is calculated as average of two central slopes (tan(angle)).\cr
#' @return Matrix of estimates and confidence intervals for intercept and slope. No standard errors provided by this algorithm.  
mc.paba <- function(angM, X, Y, alpha=0.05, posCor=TRUE, calcCI=TRUE, slope.measure=c("radian","tangent")) 
{
	## Check validity of parameters
    stopifnot(nrow(angM)==ncol(angM))
    stopifnot(length(X)==length(Y))
    stopifnot(length(X)==nrow(angM))
    stopifnot(!is.na(posCor))
	
	## Setup
	nData <- nrow(angM)
	## Valid mini slopes
	nAllItems <- sum(!is.na(angM))
	## Borderline cases at -Pi/4 or Pi/4
	nNeg <- sum(angM<=(-pi/4),na.rm=T) # nNeg = slopes <= -1 
	nNeg2 <- sum(angM<(-pi/4),na.rm=T) # nNeg2 = slopes < -1
	nPos <- sum(angM>=(pi/4),na.rm=T)
	nPos2 <- sum(angM>(pi/4),na.rm=T)
	
	##
	## Slope
	##
	## offset depending on correlation
	## for BaPa offset all with slope < -1 (slope > 1 if neg)
	## because we have also added the slopes = +/- 1 to the array
	## we have to add to the offset also half of the number of the slopes +/-.
	## the function returns twice of this offset so
	## offset = 2*(nNeg2 + (nNeg-nNeg2)/2) = nNeg + nNeg2
	if(posCor) 
        nOffset2 <- nNeg+nNeg2
	else 
        nOffset2 <- -nPos-nPos2
	## Two times index of median
	nValIndex2 <- nAllItems + nOffset2	
	## Sort and select median (NAs are dropped)
	half <- (nValIndex2+1L)%/%2L
    
    ## Extreme case half==0 can happen and results in invalid index => simply shift by 1
    if(half==0) half <- 1

    ## Shifted median should always be valid data index
    stopifnot(half<=nAllItems & half>0)

    ## Extreme case half==nAllItems can happen, for even nAllItems equation for b in PaBa paper
    ## then undefined since second index out of bounds => simply calculate b as in odd nAllItems case
	if((nValIndex2%%2L == 1L) || (half==nAllItems)) 
    {
		if(calcCI) 
            sortedM <- sort(angM) # complete sorting
		else 
            sortedM <- sort(angM,partial=half) # only partial sorting
		slope <- sortedM[half]
	} 
	else 
    {
		if(calcCI) 
            sortedM <- sort(angM) # complete sorting
		else 
            sortedM <- sort(angM,partial=half+0L:1L) # only partial sorting
        
        
        if(slope.measure=="radian") slope <- mean(sortedM[half+0L:1L])
        else slope <- atan(mean(tan(sortedM[half+0L:1L]))) # slope.measure = tangent
	}
	mcres.slope <- tan(slope)
	
    ##
	## Confidence intervals for slope
	##
	if(calcCI) 
    {
		dConf <- qnorm(1-alpha/2)*sqrt(nData*(nData-1)*(2*nData+5)/18)
		## Lower CI
		## NOTE: Difference to original Passing & Bablok paper, there (nAllItems-dConf)/2 is rounded.
		##       No averaging in case of even indices is performed. Here we calculate the quantile in analogy
		##       to the median, which seems more accurate.
		nInd <- round(nAllItems-dConf+nOffset2)
		if(posCor) 
            LowestIdx <-  2*(nNeg-nNeg2)+1
		else 
            LowestIdx <-  2*(nPos-nPos2)+1
		if(nInd < LowestIdx) 
            mcres.slopeL = -Inf
		else 
        {
			half <- (nInd+1L)%/%2L
			if(nInd%%2L == 1L) 
                slope <- sortedM[half]
			else {
                if(slope.measure=="radian") slope <- mean(sortedM[half+0L:1L])
                else slope <- atan(mean(tan(sortedM[half+0L:1L]))) # slope.measure = tangent
            }
                
			mcres.slopeL = tan(slope)
		}
		## Upper CI	
		nInd  <- round(nAllItems+dConf+nOffset2)
		if(nInd > 2*nAllItems-1) 
            mcres.slopeU <- Inf
		else 
        {
			half <- (nInd+1L)%/%2L
			if(nInd%%2L == 1L) 
                slope <- sortedM[half]
            else {
                if(slope.measure=="radian") slope <- mean(sortedM[half+0L:1L])
                else slope <- atan(mean(tan(sortedM[half+0L:1L]))) # slope.measure = tangent
            }
            mcres.slopeU = tan(slope)
		}
	}
	else 
    {
		## No theoretical CIs computed, use resampling to get CIs
		mcres.slopeL <- as.numeric(NA)
		mcres.slopeU <- as.numeric(NA)
	}
	
	##
	## Intercept
	##
	mcres.intercept <- median(calcDiff(Y,mcres.slope*X))
	if(calcCI) 
    {
		if(mcres.slopeL==-Inf) 
            mcres.interceptL <- -Inf
		else 
            mcres.interceptL <- median(calcDiff(Y,mcres.slopeU*X))
        
		if(mcres.slopeU==Inf) 
            mcres.interceptU <- Inf
		else 
            mcres.interceptU <- median(calcDiff(Y,mcres.slopeL*X))
	}
	else 
    {
		mcres.interceptL <- as.numeric(NA)
		mcres.interceptU <- as.numeric(NA)
	}
	
	## Prepare result matrix
	rmat <- matrix(nrow=2,ncol=4)
	rownames(rmat) <- c("Intercept","Slope")
	colnames(rmat) <- c("EST","SE","LCI","UCI")
	rmat[,1] <- c(mcres.intercept,mcres.slope)
	rmat[,2] <- NA
	rmat["Intercept","LCI"] <- mcres.interceptL
	rmat["Intercept","UCI"] <- mcres.interceptU
	rmat["Slope","LCI"] <- mcres.slopeL
	rmat["Slope","UCI"] <- mcres.slopeU
	return(rmat)
}
	
