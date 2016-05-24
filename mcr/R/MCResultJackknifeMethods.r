###############################################################################
##
## MCResultJackknifeMethods.R
##
## Definition of methods for class MCResultJackknife
## Class of mcreg result objects that contain Jackknife based results.
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

#' MCResultJackknife Object Constructor with Matrix in Wide Format as Input
#'
#' @param wdata measurement data in matrix format. First column reference method (x), second column test method (y).
#' @param para regression parameters in matrix form. Rows: Intercept, Slope. Cols: EST, SE, LCI, UCI.
#' @param sample.names names of individual data points, e.g. barcodes of measured samples.
#' @param method.names names of reference and test method.
#' @param regmeth name of statistical method used for regression.
#' @param glob.coef global coefficients
#' @param cimeth name of statistical method used for computing confidence intervals.
#' @param B0jack jackknife intercepts
#' @param B1jack jeckknife slopes
#' @param error.ratio ratio between standard deviation of reference and test method.
#' @param alpha numeric value specifying the 100(1-\code{alpha})\% confidence level of confidence intervals (Default is 0.05).
#' @param weight numeric vector specifying the weights used for each point
#' @return MCResult object containing regression results.
newMCResultJackknife <- function(   wdata, para, sample.names=NULL, method.names=NULL, regmeth="Unknown", glob.coef,
                                    cimeth="unknown",B0jack,B1jack,error.ratio=error.ratio,alpha=0.05,
                                    weight=rep(1,nrow(wdata))) 
{
    ## Check validity of parameters
    stopifnot(is.matrix(wdata))
    stopifnot(ncol(wdata)==2)
    stopifnot(is.matrix(para))
    stopifnot(all(dim(para)==c(2,4)))
    stopifnot(is.character(regmeth))
    stopifnot(is.element(regmeth,c("LinReg","WLinReg","Deming","PaBa","WDeming", "PaBaLarge")))
    stopifnot(is.character(cimeth))
    stopifnot(cimeth=="jackknife")
    stopifnot(is.numeric(alpha))
    stopifnot(alpha<=1 & alpha>=0)
    stopifnot(is.numeric(glob.coef))
    stopifnot(length(glob.coef)==2)
    stopifnot(is.numeric(B0jack))
    stopifnot(is.numeric(B1jack))
    stopifnot(length(B0jack)==length(B1jack))
    stopifnot(is.numeric(error.ratio))
    stopifnot(error.ratio>0)
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
        mnames <- paste("Method", 1:2, sep="")
    else 
    {
        stopifnot(length(method.names)==2)
        mnames <- method.names
    }
    ## Build object
    data <- data.frame(sid=snames, x=wdata[,1], y=wdata[,2])
    rownames(para) <- c("Intercept","Slope")
    colnames(para) <- c("EST","SE","LCI","UCI")
    names(weight) <- snames
    
    new(Class="MCResultJackknife", data=data, para=para, mnames=mnames, regmeth, cimeth=cimeth, alpha=alpha,
        glob.coef=glob.coef, B0jack=B0jack, B1jack=B1jack, error.ratio=error.ratio, weight=weight)
}


###############################################################################
## Methods
###############################################################################

#' Initialize Method for 'MCResultJackknife' Objects.
#' 
#' Method initializes newly created objects of class 'MCResultAnalytical'.
#' 
#' @param .Object object to be initialized
#' @param data empty data.frame
#' @param para empty coefficient matrix
#' @param mnames empty method names vector
#' @param regmeth string specifying the regression-method 
#' @param cimeth string specifying the confidence interval method
#' @param error.ratio for deming regression 
#' @param alpha value specifying the 100(1-alpha)% confidence-level
#' @param glob.coef global coefficients
#' @param B0jack jackknife intercepts
#' @param B1jack jackknife slopes
#' @param weight 1 for each data point

MCResultJackknife.initialize <- function(   .Object, data=data.frame(X=NA, Y=NA), para=matrix(NA, ncol=4, nrow=2),
                                            mnames=c("unknown","unknown"),regmeth="unknown",cimeth="jackknife",alpha=0.05,
                                            glob.coef=c(0,0), B0jack=0,B1jack=0,error.ratio=0, weight=1) 
{    
    .Object@data <- data
    .Object@para <- para
    .Object@mnames <- mnames
    .Object@regmeth <- regmeth
    .Object@alpha <- alpha
    .Object@cimeth <-cimeth
    .Object@glob.coef <- glob.coef
    .Object@error.ratio <- error.ratio
    .Object@B0jack <- B0jack
    .Object@B1jack <- B1jack
    .Object@weight<-weight    
    
    return(.Object)
}


#' Caluculate Response 
#' 
#' Calculate predicted values for given values of the reference-method.
#' 
#' @param .Object object of class 'MCResultJackknife'
#' @param x.levels numeric vector specifying values of the reference method for which prediction should be made
#' @param alpha significance level for confidence intervals

MCResultJackknife.calcResponse <- function(.Object, x.levels, alpha=0.05){
    
    stopifnot(is.numeric(alpha))
    stopifnot(alpha > 0 & alpha < 1)
    stopifnot(length(alpha) > 0)
    stopifnot(is.numeric(x.levels))
    stopifnot(!is.na(x.levels))
    stopifnot(length(x.levels) > 0)
    
    npoints<-length(.Object@data[,"x"])
    
    b0 <- .Object@glob.coef[1]
    b1 <- .Object@glob.coef[2]
    
    B0jack <- .Object@B0jack
    B1jack<- .Object@B1jack
    
    ylevels <- matrix(nrow=npoints,ncol=length(x.levels))   # ylevels[k,j]= b0(k)+b1(k)*x.level[j]
    colnames(ylevels) <- paste("X",1:length(x.levels),sep="")
    
    for(i in seq(along = x.levels)) 
        ylevels[,i] <- B0jack+B1jack*x.levels[i]
    
    mresp <- matrix(ncol=5,nrow=length(x.levels))
    colnames(mresp) <- c("X","Y","Y.SE","Y.LCI","Y.UCI")
    
    if(length(x.levels) == 0) 
    {
        return(mresp) # return empty matrix if no x values
    } 
    else 
    {
        mresp[,"X"] <- x.levels
        mresp[,"Y"] <- rep(b0, length(x.levels))+rep(b1, length(x.levels))*x.levels  # mresp[j,"Y"]=glob.b0+glob.b1*x.level[j]
        Y <- matrix(nrow=npoints, ncol=length(x.levels))
        for(i in 1:npoints) 
            Y[i,] <- mresp[,"Y"]
        delta.Y <- npoints*Y-(npoints-1)*ylevels
        mresp[,"Y.SE"] <- apply(delta.Y, 2, sd, na.rm=TRUE)/sqrt(npoints)
        mresp[mresp[,"Y.SE"]==0,"Y.SE"] <- NA
        mresp[,"Y.LCI"] <- mresp[,"Y"]-qt(1-alpha/2, npoints-2)*mresp[,"Y.SE"]
        mresp[,"Y.UCI"] <- mresp[,"Y"]+qt(1-alpha/2, npoints-2)*mresp[,"Y.SE"]
        return(mresp)
    }
}


#' Get-Method for Jackknife-Slope Value.
#' 
#' Extracts the slope value from objects of class 'MCResultJackknife'.
#' 
#' @param .Object object of class 'MCResultJackknife'
#' 
#' @return (numeric) jackknife-slope

MCResultJackknife.getJackknifeSlope<-function(.Object)
{
    return(.Object@B1jack)
}

#' Get-Method for Jackknife-Intercept Value.
#' 
#' Extracts the intercept value from objects of class 'MCResultJackknife'.
#' 
#' @param .Object object of class 'MCResultJackknife'
#' 
#' @return (numeric) jackknife-intercept

MCResultJackknife.getJackknifeIntercept<-function(.Object)
{
    return(.Object@B0jack)
}

#' Relative Jackknife Influence Function
#' 
#' Calculate the value of relative jackknife function for each observation.
#'
#' @param .Object object of class "MCResultJackknife" or "MCResultResampling".
#' @return a list of the following elements:
#' \item{slope}{numeric vector containing  the values of relative jackknife function of slope.} 
#' \item{intercept}{numeric vector containing  the values of relative jackknife function of intercept.}
#' @aliases getRJIF
#' @references Efron, B. (1990)
#'            Jackknife-After-Bootstrap Standard Errors and Influence Functions.
#'            Technical Report , \bold{N 134}.
MCResultJackknife.getRJIF<-function(.Object)
{
    smeanB0 <- mean(.Object@B0jack)
    smeanB1 <- mean(.Object@B1jack)
    npoints <- length(.Object@data[,"x"])
    UB0 <- (npoints-1)*(smeanB0-.Object@B0jack) #jackknife influence function
    UB1 <- (npoints-1)*(smeanB1-.Object@B1jack)
    SUB0 <- sqrt(sum(UB0^2)/(npoints-1))
    SUB1 <- sqrt(sum(UB1^2)/(npoints-1))
    UB0star <- UB0/SUB0
    UB1star <- UB1/SUB1
    return(list(slope=UB1star, intercept=UB0star))
}

#' Plotting the Relative Jackknife Influence Function
#' 
#' The function draws reference method vs. test method as scatter plot. 
#' Observations with high influence (relative jackknife influence function is greater than 2)
#' are highlighted as red points.
#'
#' @param .Object object of class "MCResultJackknife" or "MCResultResampling"
#' @references Efron, B. (1990)
#'             Jackknife-After-Bootstrap Standard Errors and Influence Functions.
#'             Technical Report , \bold{N 134}.
#' @aliases plotwithRJIF
#' @examples
#'     #library("mcr")
#'     data(creatinine,package="mcr")
#'     x <- creatinine$serum.crea
#'     y <- creatinine$plasma.crea
#'     # Deming regression fit.
#'     # The confidence intervals for regression coefficients
#'     # are calculated with jackknife method
#'     model <- mcreg( x,y,error.ratio=1,method.reg="Deming", method.ci="jackknife",
#'                      mref.name = "serum.crea", mtest.name = "plasma.crea", na.rm=TRUE )
#'     plotwithRJIF(model)
MCResultJackknife.plotwithRJIF<-function(.Object)
{
    RJIF <- getRJIF(.Object)
    UB0star <- RJIF$intercept
    UB1star <- RJIF$slope
    XX <- .Object@data[,"x"]
    YY <- .Object@data[,"y"]
    par(mfrow=c(1,2), mar=c(8.1, 4.1, 4.1, 2.1))
    plot(.Object)
    index <- which(abs(UB1star)>2)
    points(XX[index], YY[index], pch=16, col="red")
    text(XX[index], YY[index], labels=index, adj=c(-0.5,0), col="red", cex=0.8, font=2)
    mtext(  "red color: |RJIF*|>2 for slope", side=3, line=0.3)
    plot(.Object)
    index <- which(abs(UB0star)>2)
    points(XX[index], YY[index], pch=16, col="red")
    text(XX[index], YY[index], labels=index, adj=c(-0.5,0), col="red", cex=0.8, font=2)
    mtext(  "red color: |RJIF*|>2 for intercept", side=3, line=0.3)
}

#' Jackknife Statistics
#' 
#' Calculate jackknife mean, bias and standard error.
#'
#' @param .Object object of class "MCResultJackknife" or "MCResultResampling"
#' @return table with jackknife mean, bias and standard error for intercept and slope.
MCResultJackknife.getJackknifeStatistics <- function(.Object)
{
    res <- matrix(NA, nrow=2, ncol=4)
    colnames(res) <- c("EST","Jack.Mean","Bias","Jack.SE")
    rownames(res) <- c("Intercept","Slope")
    smeanB0 <- mean(.Object@B0jack)
    smeanB1 <- mean(.Object@B1jack)
    npoints <- length(.Object@data[,"x"])
    UB0 <- (npoints-1)*(smeanB0-.Object@B0jack) #jackknife influence function
    UB1 <- (npoints-1)*(smeanB1-.Object@B1jack)
    SEjackB0 <- sqrt(sum(UB0^2)/(npoints*(npoints-1)))  # Turkey's sd ertimator for jackknife
    SEjackB1 <- sqrt(sum(UB1^2)/(npoints*(npoints-1)))
    
    res[,"EST"] <- .Object@glob.coef
    res[,"Jack.Mean"] <- c(smeanB0, smeanB1)
    res[,"Jack.SE"] <- c(SEjackB0, SEjackB1)
    res[,"Bias"] <- res[,"Jack.Mean"]-res[,"EST"]
    return(res)
}

#' Print Regression-Analysis Summary for Objects of class 'MCResultJackknife'.
#' 
#' Functions prints a summary of the regression-analysis for objects of class 'MCResultJackknife'.
#' 
#' @param .Object object of class 'MCResultJackknife'

MCResultJackknife.printSummary<-function(.Object)
{
    regmeth<-.Object@regmeth
    regtext<-""

    if(regmeth=="LinReg") 
        regtext <- "Linear Regression"
    if(regmeth=="WLinReg") 
        regtext <- "Weighted Linear Regression"
    if(regmeth=="Deming") 
        regtext <- "Deming Regression"
    if(regmeth=="WDeming") 
        regtext <- "Weighted Deming Regression"
    if(regmeth %in% c("PaBa", "PaBaLarge")) 
        regtext <- "Passing Bablok Regression"

    cat("\n\n------------------------------------------\n\n")
    cat(paste("Reference method: ",.Object@mnames[1],"\n", sep=""))
    cat(paste("Test method:     ",.Object@mnames[2],"\n", sep=""))
    cat(paste("Number of data points:", length(.Object@data[,"x"])))
    cat("\n\n------------------------------------------\n\n")
    cat("The confidence intervals are calculated with jackknife (Linnet's) method.\n")
    cat(paste("Confidence level: ",(1-.Object@alpha)*100,"%\n", sep=""))
    if(regmeth %in% c("Deming","WDeming")) 
        cat("Error ratio:", .Object@error.ratio)
    cat("\n\n------------------------------------------\n\n")
    cat(toupper(paste(regtext,"Fit:\n\n")))
    print(getCoefficients(.Object))
    cat("\n\n------------------------------------------\n\n")
    cat("JACKKNIFE SUMMARY\n\n")
    print(getJackknifeStatistics(.Object))
    return(NULL)
}


