###############################################################################
##
## MCResultAnalyticalMethods.r
##
## Definition of methods for class MCResultAnalytical
## Class of mcreg result objects that contain analytical results.
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

###############################################################################
## Constructors
###############################################################################

#' MCResultAnalytical object constructor with matrix in wide format as input.
#'
#' @param wdata Measurement data in matrix format. First column reference method (x), second column comparator method (y).
#' @param para Regression parameters in matrix form. Rows: Intercept, Slope. Cols: EST, SE, LCI, UCI.
#' @param xmean Global (weighted) mean of x-values.
#' @param sample.names Names of individual data points, e.g. barcodes of measured samples.
#' @param method.names Names of reference and comparator method.
#' @param regmeth Name of statistical method used for regression.
#' @param cimeth Name of statistical method used for computing confidence intervals.
#' @param error.ratio Ratio between standard deviation of reference and comparator method.
#' @param alpha 1 - significance level for confidence intervals.
#' @param weight numeric vector specifying the weights used for each point
#' @return MCResultAnalytical object containing regression results.
newMCResultAnalytical <- function(wdata,para,xmean,sample.names=NULL,method.names=NULL,regmeth="Unknown",cimeth="analytical",
        error.ratio=error.ratio,alpha=0.05, weight=rep(1,nrow(wdata)))
{    
    ## Check validity of parameters
    stopifnot(is.matrix(wdata))
    stopifnot(ncol(wdata)==2)
    stopifnot(is.matrix(para))
    stopifnot(all(dim(para)==c(2,4)))
    stopifnot(is.character(regmeth))
    stopifnot(is.element(regmeth,c("LinReg","WLinReg","Deming","PaBa","WDeming", "PaBaLarge")))
    stopifnot(is.character(cimeth))
    stopifnot(cimeth=="analytical")
    stopifnot(is.numeric(alpha))
    stopifnot(alpha<=1 & alpha>=0)
    stopifnot(is.numeric(xmean))
    stopifnot(is.numeric(error.ratio))
    stopifnot(error.ratio>=0)
    stopifnot(is.numeric(weight))
    stopifnot(length(weight) == nrow(wdata)) 
   
    ## Sample names
    if(is.null(sample.names))
        snames <- paste("S",1:nrow(wdata),sep="")
    else
    {
        stopifnot(length(sample.names)==nrow(wdata))
        snames <- sample.names
    }
    ## Method names
    if(is.null(method.names))
        mnames <- paste("Method",1:2,sep="")
    else 
    {
        stopifnot(length(method.names)==2)
        mnames <- method.names
    }
    ## Build object
    data <- data.frame(sid=snames,x=wdata[,1],y=wdata[,2])
    rownames(para) <- c("Intercept","Slope")
    colnames(para) <- c("EST","SE","LCI","UCI")
    names(weight) <- snames
    new(Class="MCResultAnalytical",data=data,xmean=xmean,para=para,mnames=mnames,
         regmeth,cimeth=cimeth,error.ratio=error.ratio,alpha=alpha, weight=weight)
}


###############################################################################
## Methods
###############################################################################

#' Initialize Method for 'MCResultAnalytical' Objects.
#' 
#' @param .Object object to be initialized
#' @param data empty data.frame
#' @param xmean mean value
#' @param para empty coefficient matrix
#' @param mnames empty method names vector
#' @param regmeth string specifying the regression-method 
#' @param cimeth string specifying the confidence interval method
#' @param error.ratio for deming regression 
#' @param alpha value specifying the 100(1-alpha)% confidence-level
#' @param weight 1 for each data point

MCResultAnalytical.initialize <- function(.Object,data=data.frame(X=NA,Y=NA),xmean=0,para=matrix(NA,ncol=4,nrow=2),
                                          mnames=c("unknown","unknown"),regmeth="unknown",cimeth="analytical",
                                          error.ratio=0,alpha=0.05, weight=1)
{    
    .Object@data <- data
    .Object@para <- para
    .Object@xmean <- xmean
    .Object@mnames <- mnames
    .Object@regmeth <- regmeth
    .Object@alpha <- alpha
    .Object@error.ratio<-error.ratio
    .Object@cimeth<-cimeth
    .Object@weight<-weight
        
    return(.Object)
}


#' Caluculate Response 
#' 
#' Calculate predicted values for given values of the reference-method.
#' @param .Object object of class 'MCResultAnalytical'
#' @param x.levels numeric vector specifying values of the reference method for which prediction should be made
#' @param alpha significance level for confidence intervals

MCResultAnalytical.calcResponse <- function(.Object,x.levels,alpha=0.05)
{    
    stopifnot(is.numeric(x.levels))
    stopifnot(!is.na(x.levels))
    stopifnot(length(x.levels) > 0)
    stopifnot(!is.na(alpha))
    stopifnot(is.numeric(alpha))
    stopifnot(alpha > 0 & alpha < 1)
    
    mresp <- matrix(nrow=length(x.levels),ncol=5)
    colnames(mresp) <- c("X","Y","Y.SE","Y.LCI","Y.UCI")
    npoints<-length(.Object@data[,"x"])
    
    xw <- .Object@xmean
    
    if(length(x.levels) == 0)
        return(mresp) # return empty matrix if no x values
    else
    {
        rownames(mresp) <- paste("X",1:length(x.levels),sep="")
        mresp[,"X"] <- x.levels
        mresp[,"Y"] <- .Object@para["Intercept","EST"] + .Object@para["Slope","EST"]*mresp[,"X"]

        if ( !.Object@regmeth %in% c("PaBa", "PaBaLarge"))
        {
            ## Calculation based on SE estimates for intercept and slope
            mresp[,"Y.SE"] <-
                    sqrt(.Object@para["Intercept","SE"]^2 +
                         .Object@para["Slope","SE"]^2 * mresp[,"X"] *
                         (mresp[,"X"] - 2 * rep(xw,length(mresp[,"X"]))))
            mresp[,"Y.LCI"] <- mresp[,"Y"] - qt(1-alpha/2,npoints-2)*mresp[,"Y.SE"]
            mresp[,"Y.UCI"] <- mresp[,"Y"] + qt(1-alpha/2,npoints-2)*mresp[,"Y.SE"]
        }
        else
        {
            ## No analytical estimation for PaBa available
            mresp[,c("Y.SE","Y.LCI","Y.UCI")] <- NA
        }
                
        return(mresp)
    }
}

#' Print Regression-Analysis Summary for Objects of class 'MCResultAnalytical'.
#' 
#' Function prints a summary of the regression-analysis for objects of class 'MCResultAnalytical'.
#' 
#' @param .Object object of class 'MCResultAnalytical'

MCResultAnalytical.printSummary<-function(.Object)
{
    regmeth<-.Object@regmeth
    regtext<-""
    if (regmeth=="LinReg") regtext<-"Linear Regression"
    if (regmeth=="WLinReg") regtext<-"Weighted Linear Regression"
    if (regmeth=="Deming") regtext<-"Deming Regression"
    if (regmeth=="WDeming") regtext <- "Weighted Deming Regression"
    if (regmeth %in% c("PaBa", "PaBaLarge")) regtext <- "Passing-Bablok Regression"
    cat("\n\n------------------------------------------\n\n")
    cat(paste("Reference method: ",.Object@mnames[1],"\n",sep=""))
    cat(paste("Test method:     ",.Object@mnames[2],"\n",sep=""))
    cat(paste("Number of data points:", length(.Object@data[,"x"])))
    cat("\n\n------------------------------------------\n\n")
    cat("The confidence intervals are calculated with analytical method.\n")
    cat(paste("Confidence level: ",(1-.Object@alpha)*100,"%\n",sep=""))
    if (regmeth %in% c("Deming","WDeming")) cat("Error ratio:",.Object@error.ratio)
    cat("\n\n------------------------------------------\n\n")
    cat(toupper(paste(regtext,"Fit:\n\n")))
    print(getCoefficients(.Object))
    return(NULL)
}
