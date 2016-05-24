###############################################################################
##
## mcBootstrap.R
##
## Function for computing point estimations and standard  errors
## for regression coefficients with bootstrap or jackknife method.
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

#' Resampling estimation of regression parameters and standard  errors.
#' 
#' Generate jackknife or (nested-) bootstrap replicates of a statistic applied to data. 
#' Only a nonparametric ballanced design is possible. For each sample calculate 
#' point estimations and standard  errors for regression coefficients.
#' 
#' @param X Measurement values of reference method
#' @param Y Measurement values of test method
#' @param error.ratio Ratio between squared measurement errors of reference- and test method, 
#'                    necessary for Deming regression. Default 1.
#' @param method.reg Regression method.  It is possible to choose between five regression types:
#'                        \code{"LinReg"} - ordinary least square regression, 
#'         \code{"WLinReg"} - weighted ordinary least square regression,\code{"Deming"} - Deming regression, 
#'         \code{"WDeming"} - weighted Deming regression, \code{"PaBa"} - Passing-Bablok regression. 
#' @param bootstrap Bootstrap based confidence interval estimation method. 
#' @param jackknife Logical value. If TRUE - Jackknife based confidence interval estimation method.
#' @param nsamples Number of bootstrap samples.
#' @param nnested Number of nested bootstrap samples.
#' @param iter.max maximum number of iterations for weighted Deming iterative algorithm.
#' @param threshold Numerical tolerance for weighted Deming iterative algorithm convergence.
#' @param NBins number of bins used when 'reg.method="PaBaLarge"' to classify each slope in one of 'NBins' bins of constant slope angle covering the range of all slopes.
#' @param slope.measure angular measure of pairwise slopes used for exact PaBa regression (see \code{\link{mcreg}} for details).\cr   
#'          \code{"radian"} - for data sets with even sample numbers median slope is calculated as average of two central slope angles.\cr
#'          \code{"tangent"} - for data sets with even sample numbers median slope is calculated as average of two central slopes (tan(angle)).\cr
#' @return a list consisting of 
#'  \item{glob.coef}{Numeric vector of length two with global point estimations of intercept and slope.} 
#'  \item{glob.sigma}{Numeric vector of length two with global estimations of standard errors of intercept and slope.} 
#'  \item{xmean}{Global (weighted-)average of reference method values.}
#'  \item{B0jack}{Numeric vector with point estimations of intercept for jackknife samples.
#'                The i-th element contains point estimation for data set without i-th observation} 
#'  \item{B1jack}{Numeric vector with point estimations of slope for jackknife samples.
#'                The i-th element contains point estimation for data set without i-th observation} 
#'  \item{B0}{Numeric vector with point estimations of intercept for each bootstrap sample.
#'            The i-th element contains point estimation for i-th bootstrap sample.}
#'  \item{B1}{Numeric vector with point estimations of slope for each bootstrap sample. 
#'            The i-th element contains point estimation for i-th bootstrap sample.} 
#'  \item{MX}{Numeric vector with point estimations of (weighted-)average of reference method values for each bootstrap sample. 
#'            The i-th element contains point estimation for i-th bootstrap sample.}
#'  \item{sigmaB0}{Numeric vector with estimation of standard error of intercept for each bootstrap sample. 
#'                 The i-th element contains point estimation for i-th bootstrap sample.}
#'  \item{sigmaB1}{Numeric vector with estimation of standard error of slope for each bootstrap sample. 
#'                 The i-th element contains point estimation for i-th bootstrap sample.} 
#'  \item{nsamples}{Number of bootstrap samples.}
#'  \item{nnested}{Number of nested bootstrap samples.}
#'  \item{cimeth}{Method of confidence interval calculation (bootstrap).}
#'  \item{npoints}{Number of observations.}
#' @references  Efron, B., Tibshirani, R.J. (1993)
#'              \emph{An Introduction to the Bootstrap}. Chapman and Hall.
#'              Carpenter, J., Bithell, J. (2000)
#'              Bootstrap confidence intervals: when, which, what? A practical guide for medical statisticians.
#'              \emph{Stat Med}, \bold{19 (9)}, 1141--1164.
#' @author Ekaterina Manuilova \email{ekaterina.manuilova@@roche.com}, Fabian Model \email{fabian.model@@roche.com}

mc.bootstrap <- function(method.reg=c("LinReg","WLinReg","Deming","WDeming","PaBa", "PaBaLarge"), X, Y, error.ratio,
                         nsamples=1000, jackknife=TRUE, bootstrap=c("none","bootstrap", "nestedbootstrap"),
                         nnested=25, iter.max=30, threshold=0.00000001, NBins=1000000, slope.measure=c("radian","tangent")) 
{
	method.reg <- match.arg(method.reg)
	bootstrap <- match.arg(bootstrap)
    slope.measure <- match.arg(slope.measure)

	# Check validity of parameters
	stopifnot(is.numeric(X))
    stopifnot(is.numeric(Y))
    stopifnot(length(X) == length(Y))
    stopifnot(!is.na(X) & !is.na(Y))
    stopifnot(length(X) > 0)
    
    if(method.reg %in% c("Deming", "WDeming"))
    {
        stopifnot(!is.na(error.ratio))
        stopifnot(is.numeric(error.ratio))
        stopifnot(error.ratio >= 0)
        stopifnot(length(error.ratio) > 0) 
        
        if( method.reg == "WDeming")
        {
            stopifnot(!is.na(threshold))
            stopifnot(is.numeric(threshold))
            stopifnot(threshold >= 0)
            stopifnot(length(threshold) > 0)
            stopifnot(!is.na(iter.max))
            stopifnot(is.numeric(iter.max))
            stopifnot(length(iter.max) > 0)
            stopifnot(round(iter.max) == iter.max)
            stopifnot(iter.max > 0)
        }
    }
    if( bootstrap %in% c("bootstrap", "nestedbootstrap") )
    {
        stopifnot(!is.na(nsamples))
        stopifnot(is.numeric(nsamples))
        stopifnot(round(nsamples)==nsamples)
        stopifnot(nsamples>0)
        stopifnot(length(nsamples) > 0)
        
        if( bootstrap == "nestedbootstrap")
        {
            stopifnot(!is.na(nnested))
            stopifnot(is.numeric(nnested))
            stopifnot(round(nnested)==nnested)
            stopifnot(nnested>0)
            stopifnot(length(nnested) > 0)
        }
    }
    stopifnot(!is.na(jackknife))
    stopifnot(is.logical(jackknife))
    stopifnot(length(jackknife) > 0)

	if(method.reg %in% c("WLinReg","WDeming") & (min(X)<=0 | min(Y)<=0)) 
        stop("Weighted regression for non-positive values is not available.")
	if(method.reg %in% c("PaBa", "PaBaLarge") & (min(X)<0 | min(Y)<0)) 
        stop("Passing-Bablok regression for negative values is not available.")	
		
	# Number of data points
	npoints <- length(X)

	# Choice of arguments for regfun
	if(method.reg == "LinReg")
    {
		callfun.reg <- function(idx)
        {
            if(sd(X[idx])<=0)
            {
                ## All identical points, no estimation possible
                warning("Resampling LinReg regression: all x coordinates in subsample identical - regression coefficients undetermined!")
                return(list(b0=as.numeric(NA),b1=as.numeric(NA),se.b0=as.numeric(NA),se.b1=as.numeric(NA),xw=as.numeric(NA)))
            }
            return(mc.linreg(X[idx],Y[idx]))
        }
    }
	else if(method.reg == "WLinReg")
    {
		callfun.reg <- function(idx)
        {
            if(sd(X[idx])<=0)
            {
                ## All identical points, no estimation possible
                warning("Resampling WLinReg regression: all x coordinates in subsample identical - regression coefficients undetermined!")
                return(list(b0=as.numeric(NA),b1=as.numeric(NA),se.b0=as.numeric(NA),se.b1=as.numeric(NA),xw=as.numeric(NA)))
            }
            return(mc.wlinreg(X[idx],Y[idx]))
        }
    }
	else if(method.reg == "Deming")
    {
		callfun.reg <- function(idx)
        {
            if(sd(X[idx])<=0)
            {
                ## All identical points, no estimation possible
                warning("Resampling Deming regression: all x coordinates in subsample identical - regression coefficients undetermined!")
                return(list(b0=as.numeric(NA),b1=as.numeric(NA),se.b0=as.numeric(NA),se.b1=as.numeric(NA),xw=as.numeric(NA)))
            }
            return(mc.deming(X[idx],Y[idx],error.ratio=error.ratio))
        }
    }
	else if(method.reg == "WDeming")
    {
		callfun.reg <- function(idx)
        {
            if(sd(X[idx])<=0)
            {
                ## All identical points, no estimation possible
                warning("Resampling WDeming regression: all x coordinates in subsample identical - regression coefficients undetermined!")
                return(list(b0=as.numeric(NA),b1=as.numeric(NA),iter=0,xw=as.numeric(NA)))
            }
        	return(mc.wdemingConstCV(X[idx],Y[idx],error.ratio=error.ratio,iter.max=iter.max,threshold=threshold))
        }
    }
    else if(method.reg == "PaBaLarge")
    {
        callfun.reg <- function(idx)
        {
            if((sd(X[idx])<=0)|(sd(Y[idx])<=0))
            {
                ## All identical points on one axis, no estimation possible
                warning("Resampling PaBaLarge regression: x or y coordinates of all data points in subsample identical - regression coefficients undetermined!")
                return(list(b0=as.numeric(NA),b1=as.numeric(NA),xw=as.numeric(NA)))
            }
            paba.posCor <- cor(X[idx], Y[idx], method="kendall") >= 0
            tmpRes <- mc.paba.LargeData(X[idx], Y[idx], NBins=NBins, posCor=paba.posCor)
            return( list( b0=tmpRes["Intercept", "EST"], b1=tmpRes["Slope", "EST"], xw=as.numeric(NA) ) )
        }
    }
	else if(method.reg == "PaBa") 
    {
        ## For slope matrix determine if slope 1 or -1 is expected based on full data set
        paba.posCor.global <- cor(X,Y,method="kendall") >= 0
		## Compute slope matrix once for all further computations,
        ## global estimation of posCor means global +/-Inf assignment - needs to be corrected if correlation changes due to resampling
        paba.angM.global <- mc.calcAngleMat(X,Y,posCor=paba.posCor.global)
        ## Regression function
        callfun.reg<- function(idx)
        {
            if((sd(X[idx])<=0)|(sd(Y[idx])<=0))
            {
                ## All identical points on one axis, no estimation possible
                warning("Resampling PaBa regression: x or y coordinates of all data points in subsample identical - regression coefficients undetermined!")
                return(list(b0=as.numeric(NA),b1=as.numeric(NA),xw=as.numeric(NA)))
            }

            ## For actual regression determine if slope 1 or -1 is expected based on resampled data set!
            paba.posCor <- cor(X[idx],Y[idx],method="kendall") >= 0
            ## Run Passing-Bablok
		  	if(paba.posCor.global==paba.posCor)
            {
                ## Correlation sign did not change by resampling - use resampled angle matrix
                mc.res <- mc.paba(paba.angM.global[idx,idx],X[idx],Y[idx],posCor=paba.posCor,calcCI=FALSE,slope.measure=slope.measure)
            } else {
                ## Correlation sign did change by resampling - re-calculate angle matrix
                paba.angM.rec <- mc.calcAngleMat(X[idx],Y[idx],posCor=paba.posCor)
                mc.res <- mc.paba(paba.angM.rec,X[idx],Y[idx],posCor=paba.posCor,calcCI=FALSE,slope.measure=slope.measure)
                rm(paba.angM.rec)
            }
			return(list(b0=mc.res["Intercept","EST"],b1=mc.res["Slope","EST"],xw=as.numeric(NA)))
        }
    }
	else stop("Unknown regression function!")

	# Calculation of regression parameters for whole data set (global coefficients)
	d <- callfun.reg(1:npoints)
	glob.b0 <- d$b0
	glob.b1 <- d$b1
	xmean <- d$xw
	if (method.reg %in% c("PaBa", "PaBaLarge"))
    { 
        weight <- rep(1, npoints) 
    } 
    else 
    {  
        weight <- d$weight
    } 
	if(method.reg %in% c("WDeming","PaBa", "PaBaLarge")) 
    {
		## Analytical standard errors not available for weighted Deming and Passing-Bablok
		glob.seb0 <- as.numeric(NA)
		glob.seb1 <- as.numeric(NA)
	} 
    else 
    {
		glob.seb0 <- d$se.b0
		glob.seb1 <- d$se.b1
	}
    rm(d)

  # ----------------------------------------------------------------------------

	## Calculation of jackknife coefficients
	
	if(jackknife==TRUE) 
    {	
	    B0jack <- c()
        B1jack <- c()
                    
    	for(i in 1:npoints) 
        {
            ## Draw bootstrap (jackknife)sample
    		index <- (1:npoints)[-i]
    		## Regression
    		d<-callfun.reg(index)
    		B0jack[i] <- d$b0
    		B1jack[i] <- d$b1
    	} # end jackknife
    
   	    if(method.reg == "WDeming") 
        {		
            ## Use Linnet's algorithm to estimate parameter SEs	
            jackLinnetB0 <- mc.calcLinnetCI(B0jack,glob.b0,0.05)
            jackLinnetB1 <- mc.calcLinnetCI(B1jack,glob.b1,0.05)
            glob.seb0<-jackLinnetB0$se
            glob.seb1<-jackLinnetB1$se
        }
	
	    if(bootstrap=="none") 
        {
	        return(list(glob.coef=c(glob.b0,glob.b1),glob.sigma=c(glob.seb0,glob.seb1),xmean=xmean,
        	            B0jack=B0jack, B1jack=B1jack,npoints=npoints,cimeth="jackknife", weight=weight))
        }
    } # end jackknife==TRUE
	
	#-----------------------------------------------------------------------------
	
	if(bootstrap %in% c("bootstrap", "nestedbootstrap")) 
    {
        ##
	    ## Confidence intervals with Bootstrap or nested Bootstrap
       	##
        
        ## Vectors to save parameter estimations of b0 and b1 for each bootstrap sample
    	B0 <- B1 <- c()
    	sigmaB0 <- sigmaB1 <- c()
    	MX <- c()
    
    	## Bootstrap
    	for(i in 1:nsamples) 
        {
            ## Draw bootstrap sample
    		index <- sample(1:npoints,size=npoints,replace=TRUE)
    
    		## Regression on bootstrap sample
    		d <- callfun.reg(index)
    		## Parameter point estimates
    		B0[i] <- d$b0
    		B1[i] <- d$b1
    			
    		if(method.reg != "PaBa") 
              MX[i] <- d$xw
            else
                MX[i] <- mean(X[index])
    			
    		## Parameter SE estimation
    		if(bootstrap == "nestedbootstrap")      # Use nested bootstrap, works for all regression methods
            { 
    			B0n <- c()
    			B1n <- c()
    
    			## Bootstrap
    			for(ii in 1:nnested) 
                {
    			    ## Draw nested bootstrap sample
    				nindex <- sample(1:npoints,size=npoints,replace=TRUE)
    				## Regression on nested bootstrap sample
    				d <- callfun.reg(index[nindex])
    					B0n[ii] <- d$b0
    					B1n[ii] <- d$b1
    			} #end nest
    
                ## Estimate SD from nested bootstrap results
                sigmaB0[i] <- sd(B0n,na.rm=T)
                sigmaB1[i] <- sd(B1n,na.rm=T)
            }
            else        # Use analytical SD estimates
            { 
                if(method.reg %in% c("LinReg","WLinReg","Deming")) 
                {
                    sigmaB0[i] <- d$se.b0
                    sigmaB1[i] <- d$se.b1
                }
                else if(method.reg == "WDeming") 
                {
                    ## Use global SE estimates according to Linnet's method
                    sigmaB0[i] <- glob.seb0
                    sigmaB1[i] <- glob.seb1
                }
                else 
                    sigmaB0[i] <- sigmaB1[i] <- as.numeric(NA) # No analytical SE estimates available
            }
        } # end of bootstrap loop

        ## Check whether there were any invalid regressions due to drawing all identical samples
        na.mask <- !(is.na(B0)|is.na(B1))
        ## If more than 5% of bootstrap runs had invalid regressions => Error
        if(sum(!na.mask)/length(na.mask)>0.05) stop("There were too many resamples with undetermined regression coefficients!")
        
    	if (jackknife == FALSE) 
        {
            ## Return non-NA bootstrap results
            return(list(glob.coef=c(glob.b0,glob.b1), glob.sigma=c(glob.seb0, glob.seb1), xmean=xmean,
    	                B0=B0[na.mask], B1=B1[na.mask], MX=MX[na.mask], sigmaB0=sigmaB0[na.mask], sigmaB1=sigmaB1[na.mask],
    	                nsamples=sum(na.mask), nnested=nnested, npoints=npoints, cimeth=bootstrap, weight=weight))
        } 
        else 
        {
            if(any(is.na(B0jack)|is.na(B1jack))) stop("Too many identical data points!")
            ## Return non-NA bootstrap results
            return(list(glob.coef=c(glob.b0, glob.b1), glob.sigma=c(glob.seb0, glob.seb1), xmean=xmean,
                        B0jack=B0jack, B1jack=B1jack, B0=B0[na.mask], B1=B1[na.mask], MX=MX[na.mask], sigmaB0=sigmaB0[na.mask], sigmaB1=sigmaB1[na.mask],
              	        nsamples=sum(na.mask), nnested=nnested, npoints=npoints, cimeth=bootstrap, weight=weight))
        }
    }  # end of bootstrap
} # end of mc.bootstrap

