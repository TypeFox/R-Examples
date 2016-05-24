###############################################################################
##
## MCResultMethods.r
##
## Definition of methods for class MCResult
## Base class of mcreg result objects.
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

#' MCResult Object Constructor with Matrix in Wide Format as Input
#'
#' @param wdata measurement data in matrix format. First column reference method (x), second column test method (y).
#' @param para regression parameters in matrix form. Rows: Intercept, Slope. Cols: EST, SE, LCI, UCI.
#' @param sample.names names of individual data points, e.g. barcodes of measured samples.
#' @param method.names names of reference and test method.
#' @param regmeth name of statistical method used for regression.
#' @param cimeth name of statistical method used for computing confidence intervals.
#' @param error.ratio ratio between standard deviation of reference and test method.
#' @param alpha numeric value specifying the 100(1-\code{alpha})\% confidence level of confidence intervals (Default is 0.05).
#' @param weight numeric vector specifying the weights used for each point
#' @return MCResult object containing regression results.
newMCResult <- function(wdata, para, sample.names=NULL, method.names=NULL, regmeth="Unknown", 
                        cimeth, error.ratio, alpha=0.05, weight=rep(1,nrow(wdata))) 
{
	## Check validity of parameters
  	stopifnot(is.matrix(wdata))
	  stopifnot(ncol(wdata)==2)
    stopifnot(is.matrix(para))
    stopifnot(all(dim(para)==c(2,4)))
    stopifnot(is.character(regmeth))
    stopifnot(is.element(regmeth,c("LinReg","WLinReg","Deming","BaPa","WDeming", "PaBaLarge")))
    stopifnot(!is.na(alpha))
    stopifnot(is.numeric(alpha))
    stopifnot(length(alpha) > 0)
    stopifnot(alpha > 0 & alpha < 1)
    stopifnot(is.element(cimeth,c("analytical","jackknife","bootstrap","nestedbootstrap")))
    stopifnot(!is.na(error.ratio))
    stopifnot(is.numeric(error.ratio))
    stopifnot(length(error.ratio) > 0)
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
	new(Class="MCResult",data=data,para=para,mnames=mnames,regmeth=regmeth,cimeth=cimeth,error.ratio=error.ratio,alpha=alpha,weight=weight)
}


###############################################################################
## Methods
###############################################################################

#' MCResult Object Initialization
#'
#' @param .Object object of class "MCResult"
#' @param data measurement data in matrix format. First column reference method (x), second column test method (y).
#' @param para regression parameters in matrix form. Rows: Intercept, Slope. Cols: EST, SE, LCI, UCI.
#' @param mnames names of reference and test method.
#' @param regmeth name of statistical method used for regression.
#' @param cimeth name of statistical method used for computing confidence intervals.
#' @param error.ratio ratio between standard deviation of reference and test method.
#' @param alpha numeric value specifying the 100(1-\code{alpha})\% confidence level of confidence intervals (Default is 0.05).
#' @param weight weights to be used for observations
#' @return MCResult object with initialized parameter.
MCResult.initialize <- function(.Object, data=data.frame(X=NA,Y=NA), para=matrix(NA,ncol=4,nrow=2),
		                        mnames=c("unknown","unknown"), regmeth="unknown", cimeth="unknown", 
                                error.ratio=0,alpha=0.05, weight=1) 
{	
	.Object@data <- data
	.Object@para <- para
	.Object@cimeth <- cimeth
	.Object@mnames <- mnames
	.Object@regmeth <- regmeth
	.Object@error.ratio<-error.ratio
	.Object@alpha <- alpha
	.Object@weight <- weight
	return(.Object)
}

#' Get Regression Coefficients
#'
#' @param .Object object of class "MCResult".
#' @return Regression parameters in matrix form. Rows: Intercept, Slope. Cols: EST, SE, LCI, UCI.
#' @aliases getCoefficients

MCResult.getCoefficients <- function(.Object)
{
	return(.Object@para)
}


#' Get Error Ratio
#'
#' @param .Object Object of class "MCResult"
#' @return Error ratio. Only relevant for Deming type regressions.
#' @aliases getErrorRatio
MCResult.getErrorRatio <- function(.Object)
{
    return(.Object@error.ratio)
}

#' Get Weights of Data Points
#'
#' @param .Object Object of class "MCResult"
#' @return Weights of data points. 
#' @aliases getWeights
MCResult.getWeights <- function(.Object)
{
    return(.Object@weight)
}

#' Get Data
#'
#' @param .Object object of class "MCResult".
#' @return Measurement data in matrix format. First column contains reference method (X), second column contains test method (Y).
#' @aliases getData
MCResult.getData <- function(.Object)
{
	## Transform data frame to matrix, names could be changed
	rmat <- cbind(X=.Object@data$x,Y=.Object@data$y)
	rownames(rmat) <- .Object@data$sid
	return(rmat)
}

#' Get Regression Residuals
#' 
#' This function returnes residuals in x-direction (x-xhat),
#' in y-direction(y-yhat) and optimized residuals.
#' The optimized residuals correspond to distances between 
#' data points and the regression line which were optimized for
#' regression coefficients estimation. In case of Passing-Bablok
#' Regression orthogonal residuals will be returned as optimized residuals .
#' The residuals in x-direction are interesting for regression types which
#' assume errors in both variables  (deming, weighted deming, Passing-Bablok), 
#' particularily for checking of model assumptions.
#' @param .Object object of class "MCResult".
#' @return residuals as data frame.
#' @seealso \code{\link{plotResiduals}} 
#' @aliases getResiduals
MCResult.getResiduals <- function(.Object) 
{
    weight <- getWeights(.Object)
    x <- .Object@data$x
    y <- .Object@data$y
    lmb <- .Object@error.ratio
    a <- .Object@para["Intercept",1]
    b <- .Object@para["Slope",1]
    d <- y-(a+b*x)
           
    if (.Object@regmeth %in% c("LinReg","WLinReg")){
        xhat <- x
        yhat <- a + b*x
    } else {
        xhat <- x+lmb*b*d/(1+lmb*b^2)
        yhat <- y-d/(1+lmb*b^2)
    }
    xres <- x-xhat
    yres <- y-yhat
 
    if (.Object@regmeth %in% c("LinReg","WLinReg")){
        res <- yres*sqrt(weight)
    } else {
        res <- sign(yres)*sqrt((xres)^2+lmb*(yres)^2)*sqrt(weight)
    }     

    result <- data.frame(x=xres, y=yres, optimized=res)

    return(result)
}


#' Get Fitted Values.
#' 
#' This funcion computes fitted values for a 'MCResult'-object. Depending
#' on the regression method and the error ratio, a projection onto the regression
#' line is performed accordingly. For each point (x_i; y_i) i=1,...,n the projected
#' point(x_hat_i; y_hat_i) is computed.
#' 
#' @param .Object object of class "MCResult".
#' @return fitted values as data frame.
#' @seealso \code{\link{plotResiduals}} \code{\link{getResiduals}} 
#' @aliases getFitted
MCResult.getFitted <- function(.Object)
{
    x <- .Object@data$x
    y <- .Object@data$y
    lmb <- .Object@error.ratio
    a <- .Object@para["Intercept",1]
    b <- .Object@para["Slope",1]
    
    d <- y-(a+b*x)
        
    if (.Object@regmeth %in% c("LinReg","WLinReg"))
    {
        xhat <- x
        yhat <- a + b*x
    } 
    else
    { 
        xhat <- x+lmb*b*d/(1+lmb*b^2)
        yhat <- y-d/(1+lmb*b^2)
    } 
    
    return(data.frame(x_hat=xhat,y_hat=yhat))
}


#' Get Regression Method
#'
#' @param .Object object of class "MCResult".
#' @return Name of the statistical method used for the regression analysis.
#' @aliases getRegmethod
MCResult.getRegmethod <- function(.Object)
{
	return(.Object@regmeth)
}


###
### Statistics Functions
###

#' Calculate CUSUM Statistics According to Passing & Bablok (1983)
#'
#' @param .Object object of class "MCResult".
#' @return A list containing the following elements:
#' \item{nPos}{sum of positive residuals}
#' \item{nNeg}{sum of negative residuals}
#' \item{cusum}{a cumulative sum of vector with scores ri for each point, sorted increasing by distance of points to regression line.}
#' \item{max.cumsum}{Test statisics of linearity test}
#' @references Passing, H., Bablok, W. (1983) 
#'              A new biometrical procedure for testing the equality of measurements from two different analytical methods. 
#'              Application of linear regression procedures for method comparison studies in clinical chemistry, Part I. 
#'              \emph{J Clin Chem Clin Biochem}.  Nov; \bold{21(11)}:709--20.
#' @aliases calcCUSUM
MCResult.calcCUSUM <- function(.Object)
{
	
  res <- .Object@data[,"y"]-.Object@para["Intercept","EST"]-.Object@data[,"x"]*.Object@para["Slope","EST"]
	nl <- sum(res>0)
	nL <- sum(res<0)
	ri <- ifelse(res>0,sqrt(nL/nl),ifelse(res<0,-sqrt(nl/nL),0))
	Di <- (.Object@data$y+.Object@data$x/.Object@para["Slope",1]-.Object@para["Intercept",1])/sqrt(1+1/(.Object@para["Slope",1]^2))
	cusum <- cumsum(ri[order(Di)])
	max.cusum <- max(abs(cusum),na.rm=TRUE)
	return(list(nPos=nl,nNeg=nL,cusum=cusum,max.cusum=max.cusum))
}

#' Bland-Altman Plot
#' 
#' Draw different Bland-Altman plot modifications (see parameter \code{plot.type}).
#'
#' @param .Object object of class "MCResult".
#' @param xlab label for the x-axis
#' @param ylab label for the y-axis
#' @param digits number of decimal places for the difference of means and standard deviation appearing in the plot. 
#' @param plot.type integer specifying a specific Bland-Altman plot modification (default is 3).
#'  Possible choices are:
#'  1 - difference plot X vs. Y-X with null-line and mean plus confidence intervals.\cr
#'  2 - difference plot X vs. (Y-X)/X (relative differences) with null-line and mean.\cr
#'  3 - difference plot 0.5*(X+Y) vs. Y-X with null-line and mean plus confidence intervals.\cr
#'  4 - difference plot 0.5*(X+Y) vs. (Y-X)/X  (relative differences) with null-line.\cr
#'  5 - difference plot rank(X) vs. Y-X with null-line and mean plus confidence intervals.\cr
#'  6 - difference plot rank(X) vs. (Y-X)/X (relative differences) with null-line and mean.\cr
#'  7 - difference plot sqrt(X*Y) vs. Y/X with null-line and mean plus confidence intervals calculated with help of log-transformation.\cr
#'  8 - difference plot 0.5*(X+Y) vs. (Y-X) / (0.5*(X+Y)) with null-line.\cr
#' @param main plot title.
#' @param ref.line logical value. If \code{ref.line=TRUE} (default), the reference line will be drawn.
#' @param ref.line.col reference line color.
#' @param ref.line.lty reference line type.
#' @param ref.line.lwd reference line width.
#' @param bias.line.lty line type for estimated bias.
#' @param bias.line.lwd line width for estimated bias.
#' @param bias.line.col color of the line for estimated bias.
#' @param bias.text.col color of the label for estimated bias (defaults to the same as \code{bias.line.col}.)
#' @param bias.text.cex The magnification to be used for the label for estimated bias 
#' @param loa.line.lty line type for estimated limits of agreement.
#' @param loa.line.lwd line width for estimated limits of agreement.
#' @param loa.line.col color of the line for estimated limits of agreement.
#' @param loa.text.col color of the label for estimated limits of agreement (defaults to the same as \code{loa.line.col}.)
#' @param add.grid logical value. If \code{add.grid=TRUE} (Default) gridlines will be drawn.
#' @param ylim limits for the y-axis
#' @param cex numeric value specifying the magnification factor used for points
#' @param ... further graphical parameters
#' @references  Bland, J. M., Altman, D. G. (1986)
#'              Statistical methods for assessing agreement between two methods of clinical measurement. 
#'              \emph{Lancet}, \bold{i:} 307--310. 
#' @seealso \code{\link{plot.mcr}}, \code{\link{plotResiduals}}, \code{\link{plotDifference}}, \code{\link{plotBias}}, \code{\link{compareFit}}
#' @aliases plotDifference
#' @examples
#'     #library("mcr")
#'     data(creatinine,package="mcr")
#'     x <- creatinine$serum.crea
#'     y <- creatinine$plasma.crea
#' 
#'     # Deming regression fit.
#'     # The confidence intercals for regression coefficients
#'     # are calculated with analytical method
#'     model <- mcreg( x,y,error.ratio=1,method.reg="Deming", method.ci="analytical",
#'                      mref.name = "serum.crea", mtest.name = "plasma.crea", na.rm=TRUE )
#' 
#'     plotDifference( model ) # Default plot.type=3
#'     plotDifference( model, plot.type=5)
#'     plotDifference( model, plot.type=7, ref.line.lty=3, ref.line.col="green3" )

MCResult.plotDifference <- function(.Object,
                                    xlab=NULL,
                                    ylab=NULL,
                                    ref.line=TRUE,
                                    ref.line.col="black",
                                    ref.line.lty=1,
                                    ref.line.lwd=1,
                                    bias.line.lty=1,
                                    bias.line.lwd=1,
                                    bias.line.col="red",
                                    bias.text.col=NULL,
                                    bias.text.cex=0.8,
                                    loa.line.lty=2,
                                    loa.line.lwd=1,
                                    loa.line.col="red",
                                    loa.text.col=NULL, 
                                    plot.type=3,
                                    main=NULL,
                                    cex=0.8, 
                                    digits=2, 
                                    add.grid= TRUE,
                                    ylim=NULL,
                                    ...)
{
	stopifnot(is.element(plot.type,1:8))
	stopifnot(is.numeric(digits))
	stopifnot(digits>0)
	stopifnot(digits==round(digits))
	
	if(is.null(main)) main <- "Difference plot" 
	
	if(plot.type %in% c(1,2))
    {
		X <- .Object@data[,"x"]
		x.lab <- .Object@mnames[1]
	}
		
	if(plot.type %in% c(3,4,8))
    {
		X <- (.Object@data[,"x"]+.Object@data[,"y"])/2
		x.lab <- paste("(",.Object@mnames[1],"+",.Object@mnames[2],")/2")
	}
	
	if(plot.type == 7)
    {
        stopifnot(all(.Object@data[,"x"]*.Object@data[,"y"]>0))
        X <- sqrt(.Object@data[,"x"]*.Object@data[,"y"])
        text1 <-.Object@mnames[1]
        text2 <- .Object@mnames[2]
        text1 <- gsub('[[:punct:]]',".",text1)
        text1 <- gsub('[[:space:]]',".",text1)
        text2 <- gsub('[[:punct:]]',".",text2)
        text2 <- gsub('[[:space:]]',".",text2)
        eval(parse(text=paste("x.lab<-expression(sqrt(paste(",text1,",' * ',",text2,")))",sep="")))
	}
	                                                  
	if(plot.type %in% c(5,6))
    {                     
		X <- rank(.Object@data[,"x"])
		x.lab <- paste("Rank of",.Object@mnames[1])
	}
	
	if(plot.type %in% c(1,3,5))
    {
		Y <- .Object@data[,"y"]-.Object@data[,"x"]
		y.lab <- paste(.Object@mnames[2],"-",.Object@mnames[1])
	}
	
	if(plot.type %in% c(2,4,6))
    {
        stopifnot(all(.Object@data[,"x"]!= 0))
        Y <- (.Object@data[,"y"]-.Object@data[,"x"])/.Object@data[,"x"]
        y.lab<-paste("(",.Object@mnames[2],"-",.Object@mnames[1],")/",.Object@mnames[1])
	}
	
	if(plot.type == 7)
    {
        stopifnot(all(.Object@data[,"x"]!= 0))
        Y <- .Object@data[,"y"]/.Object@data[,"x"]
        y.lab <- paste(text2,"/",text1)
	}
	
    if(plot.type == 8)
    {
        stopifnot(all((.Object@data[,"x"]+.Object@data[,"y"])/2 != 0))
        Y <- (.Object@data[,"y"]-.Object@data[,"x"])/((.Object@data[,"y"]+.Object@data[,"x"])/2 )
        text1 <- .Object@mnames[1]
        text2 <- .Object@mnames[2]
        y.lab <- paste("( ",text2," - ",text1," )"," / mean( ",text1," , ",text2," )", sep="" )
    }
	
	if(plot.type %in% c(1,2,3,4,5,6,8))
    {
        MEAN <- mean(Y,na.rm=TRUE)
        LOW <- MEAN - sd(Y)*2
        UP <- MEAN + sd(Y)*2
        y.low <- min(min(Y,na.rm=TRUE),LOW-0.05*2*sd(Y),na.rm=TRUE)
        y.up <- max(max(Y,na.rm=TRUE),UP+0.05*2*sd(Y),na.rm=TRUE)
        YLim <- c(y.low,y.up)   ## ylim replaced by YLim
	}
	
    if(plot.type == 7)
    {
        XX <- .Object@data[,"x"]
        YY <- .Object@data[,"y"]
        logmean <- mean(log(YY/XX),na.rm=TRUE)
        logLOW <- logmean - 2*sd(log(YY/XX),na.rm=TRUE)
        logUP  <- logmean + 2*sd(log(YY/XX),na.rm=TRUE)
        MEAN <- exp(logmean)
        LOW <- exp(logLOW)
        UP <- exp(logUP)
        y.low <- min(min(Y,na.rm=TRUE),LOW-0.05*2*sd(log(YY/XX),na.rm=TRUE),na.rm=TRUE)
        y.up <- max(max(Y,na.rm=TRUE),UP+0.05*2*sd(log(YY/XX),na.rm=TRUE),na.rm=TRUE)
        YLim <- c(y.low,y.up)   ## ylim replaced by YLim
    }
    
    par(mar=c(5, 5, 4, 5) + 0.1)
    xlim <- c(range(X, na.rm=TRUE)[1],range(X, na.rm=TRUE)[2]+1/10*(range(X, na.rm=TRUE)[2]-range(X, na.rm=TRUE)[1]))
    
    if(is.null(xlab)) 
        xlab <- x.lab
    if(is.null(ylab)) 
        ylab <- y.lab
    
    if(!is.null(ylim))               # replace Min/Max-'ylim' by user-defined 'ylim' (YLimits)
    {
        if(all(is.numeric(ylim)))
        {
            if(all(!is.na(ylim)))
            {
                if(length(ylim) == 2)
                    YLim <- ylim
            }
        }
    }
    
    plot(X,Y,xlab=xlab,ylab=ylab,main=main,ylim=YLim,cex=cex, ...)
	
    if(add.grid) grid()
	
    if(ref.line==TRUE)
    {
        if(plot.type == 7)
            abline(h=1,
                   col=ref.line.col, 
                   lty=ref.line.lty, 
                   lwd=ref.line.lwd)
        else
            abline(h=0,
                   col=ref.line.col, 
                   lty=ref.line.lty, 
                   lwd=ref.line.lwd)
   }
  
    if(is.null(bias.text.col)) 
        bias.text.col <- bias.line.col
    if(is.null(loa.text.col)) 
        loa.text.col <- loa.line.col
  
    abline(h=MEAN, lty=bias.line.lty,
           lwd=bias.line.lwd,
           col=bias.line.col)
   
	meantext <- ifelse(plot.type %in% paste(c(1:6,8)),paste("MEAN"),paste("MEAN(log)"))
	mtext(at=MEAN+(y.up-y.low)/50,side=4, line=1,text=meantext,cex= bias.text.cex,col=bias.text.col,las=1)
	mtext(at=MEAN-(y.up-y.low)/50, side=4,line=1,text=paste(round(MEAN,digits=digits)),cex= bias.text.cex,col=bias.text.col,las=1)
	
	if(plot.type %in% c(1,3,5,7))
    {
        abline(h=UP,lty=loa.line.lty, 
               lwd=loa.line.lwd, 
               col=loa.line.col)
    	
        uptext <- ifelse(plot.type %in% paste(c(1:6,8)),paste("+ 2 SD"),paste("+ 2 SD(log)"))
        mtext(at=UP+(y.up-y.low)/50,side=4,line=1,text=uptext,cex= bias.text.cex,col=loa.text.col,las=1)
        mtext(at=UP-(y.up-y.low)/50,side=4,line=1,text=paste(round(UP,digits=digits)),cex= bias.text.cex,col=loa.text.col,las=1)
    	
        abline(h=LOW, lty=loa.line.lty, 
               lwd=loa.line.lwd, 
               col=loa.line.col)
    	
        lowtext <- ifelse(plot.type %in% paste(c(1:6,8)),paste("- 2 SD"),paste("- 2 SD(log)"))
        mtext(at=LOW+(y.up-y.low)/50,side=4,line=1,text=lowtext,cex= bias.text.cex,col=loa.text.col,las=1)
        mtext(at=LOW-(y.up-y.low)/50,side=4,line=1,text=paste(round(LOW,digits=digits)),cex= bias.text.cex,col=loa.text.col,las=1)
    }
}

#' Calculate Response with Confidence Interval.
#' 
#' Calculate Response \eqn{Intercept + Slope * Refrencemethod} with Corresponding Confidence Interval
#'
#' @param .Object object of class "MCResult".
#' @param x.levels a numeric vector with points for which response schould be calculated.
#' @param alpha numeric value specifying the 100(1-\code{alpha})\% confidence level of the confidence interval (Default is 0.05).
#' @param ... further parameters
#' @return response and corresponding confidence interval for each point in vector \code{x.levels}.
#' @seealso \code{\link{calcBias}}
#' @aliases calcResponse
#' @examples
#'     #library("mcr")
#'     data(creatinine,package="mcr")
#'     x <- creatinine$serum.crea
#'     y <- creatinine$plasma.crea
#'     # Deming regression fit.
#'     # The confidence intercals for regression coefficients
#'     # are calculated with analytical method
#'     model <- mcreg( x,y,error.ratio=1,method.reg="Deming", method.ci="analytical",
#'                      mref.name = "serum.crea", mtest.name = "plasma.crea", na.rm=TRUE )
#'     calcResponse(model, x.levels=c(1,2,3))
MCResult.calcResponse <- function(.Object,x.levels,alpha,...)
{
    cat("Called virtual basis function: MCResult.calcResponse\n")
}

#' Systematical Bias Between Reference Method and Test Method
#' 
#' Calculate systematical bias between reference and test methods 
#' at the decision point Xc as
#' \eqn{ Bias(Xc) = Intercept + (Slope-1) * Xc}
#' with corresponding confidence intervals.
#'
#' @param .Object object of class "MCResult".
#' @param type One can choose between absolute (default)  and proportional bias (\code{Bias(Xc)/Xc}).
#' @param percent logical value. If \code{percent = TRUE} the proportional bias will be calculated in percent.
#' @param x.levels a numeric vector with decision points for which bias schould be calculated.
#' @param alpha numeric value specifying the 100(1-\code{alpha})\% confidence level of the confidence interval (Default is 0.05).
#' @param ... further parameters
#' @return response and corresponding confidence interval for each decision point from x.levels.
#' @seealso \code{\link{plotBias}}
#' @aliases calcBias
#' @examples
#'     #library("mcr")
#'     data(creatinine,package="mcr")
#'     x <- creatinine$serum.crea
#'     y <- creatinine$plasma.crea
#' 
#'     # Deming regression fit.
#'     # The confidence intervals for regression coefficients
#'     # are calculated with analytical method
#'     model <- mcreg( x,y,error.ratio = 1,method.reg = "Deming", method.ci = "analytical",
#'
#'                      mref.name = "serum.crea", mtest.name = "plasma.crea", na.rm=TRUE )
#'     # Now we calculate the systematical bias
#'     # between the testmethod and the reference method
#'     # at the medical decision points 1, 2 and 3 
#'
#'     calcBias( model, x.levels = c(1,2,3))
#'     calcBias( model, x.levels = c(1,2,3), type = "proportional")
#'     calcBias( model, x.levels = c(1,2,3), type = "proportional", percent = FALSE)
MCResult.calcBias <- function(.Object, x.levels, type = c("absolute", "proportional"), percent=TRUE, alpha=0.05, ...) 
{	
    stopifnot(is.logical(percent))
    stopifnot(!is.na(x.levels))
    stopifnot(is.numeric(x.levels))
    stopifnot(length(x.levels) > 0)
    stopifnot(!is.na(alpha))
    stopifnot(is.numeric(alpha))
    stopifnot(length(alpha) > 0)
    stopifnot(alpha>=0 & alpha<=1)
    type <- match.arg(type)
	  
	## Bias calculation
	bias <- calcResponse(.Object,x.levels,alpha=alpha,...)
	bias[,c("Y","Y.LCI","Y.UCI")] <- bias[,c("Y","Y.LCI","Y.UCI")]-bias[,"X"]
	
	if ( type == "proportional" )
	{
	    stopifnot(all(x.levels != 0))
	    bias[,c("Y","Y.SE","Y.LCI","Y.UCI")] <- bias[,c("Y","Y.SE","Y.LCI","Y.UCI")]/x.levels
	    if (percent == TRUE)
      { 
             bias[,c("Y","Y.SE","Y.LCI","Y.UCI")] <- bias[,c("Y","Y.SE","Y.LCI","Y.UCI")]*100
             colnames(bias) <- c("Level","Prop.bias(%)","SE","LCI","UCI")
      } else colnames(bias) <- c("Level","Prop.bias","SE","LCI","UCI")       
  } else colnames(bias) <- c("Level","Bias","SE","LCI","UCI")
		
	## Return bias object
	rownames(bias) <- paste("X",1:length(x.levels),sep="")
	return(bias)
}

#' Scatter Plot Method X vs. Method Y
#' 
#' Plot method X (reference) vs. method Y (test) with (optional) line of identity,
#' regression line and confidence bounds for response.
#'
#' @param x object of class "MCResult".
#' @param alpha numeric value specifying the 100(1-\code{alpha})\% confidence bounds.
#' @param xn number of points (default 20) for calculation of confidence bounds.
#' @param draw.points logical value. If \code{draw.points=TRUE}, the data points will be drawn. 
#' @param xlim limits of the x-axis. If \code{xlim=NULL} the x-limits will be calculated automatically.
#' @param ylim limits of the y-axis. If \code{ylim=NULL} the y-limits will be calculated automatically.
#' @param x.lab label of x-axis. Default is the name of reference method.
#' @param y.lab label of y-axis. Default is the name of test method.
#' @param equal.axis logical value. If \code{equal.axis=TRUE} x-axis will be equal to y-axis.
#' @param add logical value. If \code{add=TRUE}, the plot will be drawn in current graphical window.
#' @param points.col Color of data points. 
#' @param points.pch Type of data points (see \code{par()}). 
#' @param points.cex Size of data points (see \code{par()}).
#' @param reg Logical value. If \code{reg=TRUE}, the regression line will be drawn.
#' @param reg.col Color of regression line.
#' @param reg.lty Type of regression line.
#' @param reg.lwd The width of regression line.
#' @param identity logical value. If \code{identity=TRUE} the identity line will be drawn.
#' @param identity.col The color of identity line.
#' @param identity.lty The type of identity line.
#' @param identity.lwd the width of identity line.
#' @param ci.area logical value. If \code{ci.area=TRUE} (default) the confidence area will be drawn.
#' @param ci.area.col the color of confidence area. 
#' @param ci.border logical value. If \code{ci.border=TRUE} the confidence limits will be drawn.
#' @param ci.border.col The color of confidence limits. 
#' @param ci.border.lty The line type of confidence limits. 
#' @param ci.border.lwd The line width of confidence limits. 
#' @param add.legend logical value. If \code{add.legend=FALSE} the plot will not have any legend.
#' @param legend.place The position of legend: "topleft","topright","bottomleft","bottomright". 
#' @param main String value. The main title of plot. If \code{main=NULL} it will include regression name.
#' @param sub String value. The subtitle of plot. If \code{sub=NULL} and \code{ci.border=TRUE} or \code{ci.area=TRUE} it will include the art of confidence bounds calculation.
#' @param add.cor Logical value. If \code{add.cor=TRUE} the correlation coefficient will be shown. 
#' @param cor.method a character string indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman", can be abbreviated.
#' @param add.grid Logical value. If \code{add.grid=TRUE} (default) the gridlines will be drawn.  
#' @param ... further graphical parameters
#' @seealso \code{\link{plotBias}}, \code{\link{plotResiduals}}, \code{\link{plotDifference}}, \code{\link{compareFit}},\code{\link{includeLegend}}
#' @aliases plot.mcr
#' @aliases plot
#' @examples
#'  library(mcr)
#'  data(creatinine,package="mcr")
#'  creatinine <- creatinine[complete.cases(creatinine),]
#'   x <- creatinine$serum.crea
#'   y <- creatinine$plasma.crea
#' 
#'   m1 <- mcreg(x,y,method.reg="Deming",  mref.name="serum.crea",
#'                                         mtest.name="plasma.crea", na.rm=TRUE)
#'   m2 <- mcreg(x,y,method.reg="WDeming", method.ci="jackknife",
#'                                         mref.name="serum.crea",
#'                                         mtest.name="plasma.crea", na.rm=TRUE)
#' 
#'   plot(m1,  xlim=c(0.5,3),ylim=c(0.5,3), add.legend=FALSE,
#'                            main="Deming vs. weighted Deming regression",
#'                            points.pch=19,ci.area=TRUE, ci.area.col=grey(0.9),
#'                            identity=FALSE, add.grid=FALSE, sub="")
#'   plot(m2, ci.area=FALSE, ci.border=TRUE, ci.border.col="red3",
#'                            reg.col="red3", add.legend=FALSE,
#'                            draw.points=FALSE,add=TRUE)
#' 
#'   includeLegend(place="topleft",models=list(m1,m2),
#'                            colors=c("darkblue","red"), design="1", digits=2)
MCResult.plot <- function(x,
                          alpha = 0.05 ,
                          xn = 20,
                          equal.axis=FALSE,
                          xlim = NULL,
                          ylim = NULL,
                          x.lab = x@mnames[1],
                          y.lab = x@mnames[2],
                          add = FALSE,
                          draw.points = TRUE,
                          points.col = "black",
                          points.pch =  1,
                          points.cex = 0.8,
                          reg = TRUE,              # regression line
                          reg.col =NULL,  
                          reg.lty =1, 
                          reg.lwd = 2,
                          identity = TRUE,         # bisecting line of an angle
                          identity.col = NULL,
                          identity.lty = 2,
                          identity.lwd = 1,
                          ci.area = TRUE,          # confidence bounds as area
                          ci.area.col = NULL, 
                          ci.border = FALSE, 
                          ci.border.col = NULL,    # confidence bounds as lines
                          ci.border.lty = 2,
                          ci.border.lwd = 1,
                          add.legend = TRUE,
                          legend.place = c("topleft","topright","bottomleft","bottomright"),
                          main = NULL,
                          sub = NULL,
                          add.cor = TRUE,
                          cor.method = c("pearson", "kendall", "spearman"),
                          add.grid= TRUE,
                          ...) 
{	
	stopifnot(is.logical(reg))
	stopifnot(is.logical(ci.area))
	stopifnot(is.logical(identity))
	stopifnot(is.logical(ci.border))
	stopifnot(is.logical(add.cor))
	stopifnot(is.logical(add.grid))
	stopifnot(is.logical(draw.points))
	stopifnot(is.logical(add.legend))
	stopifnot(is.numeric(alpha))
	stopifnot(alpha>0 & alpha<1)
	
    cor.method <- match.arg(cor.method)
    stopifnot(is.element(cor.method, c("pearson", "kendall", "spearman")))
	
    legend.place <- match.arg(legend.place)
	  stopifnot(is.element(legend.place, c("topleft","topright","bottomleft","bottomright")))
    if(length(points.col)>1) stopifnot(length(points.col)==nrow(x@data))
    if(length(points.pch)>1) stopifnot(length(points.pch)==nrow(x@data))
    
  	if(x@regmeth == "LinReg")
		titname<-"Linear Regression"
    else if(x@regmeth == "WLinReg")
		titname<-"Weighted Linear Regression"
    else if(x@regmeth == "Deming")
		titname<-"Deming Regression"
    else if(x@regmeth == "WDeming")
        titname<-"Weighted Deming Regression"
    else 
        titname <- "Passing Bablok Regression"
 
    ## Colors
    niceblue <- rgb(37/255,52/255,148/255)
    niceor <- rgb(230/255,85/255,13/255)
    niceblue.bounds <- rgb(236/255, 231/255, 242/255)
    	
    if(is.null(reg.col)) 
        reg.col <- niceblue 
    if(is.null(identity.col)) 
        identity.col <- niceor
    if(is.null(ci.area.col)) 
        ci.area.col <- niceblue.bounds
    if(is.null(ci.border.col)) 
        ci.border.col <- niceblue
	
	stopifnot(is.numeric(xn))
	stopifnot(round(xn) == xn)                     
	stopifnot(xn>1)
	
	if(is.null(xlim)) 
        rx <- range(x@data[,"x"],na.rm=TRUE)
	else
        rx <- xlim
        
	xd <- seq(rx[1],rx[2],length.out=xn)
	xd <- union(xd, rx)
	  
    ## additional points outer range for nice bounds
	delta <- abs(rx[1]-rx[2])/xn
	xd <- xd[order(xd)]
	xd.add <- c(xd[1]-delta*1:10, xd, xd[length(xd)]+delta*1:10)
		
	if(is.null(xlim)) xlim <- rx
  
	tmp.range <- range(as.vector(x@data[,c("x","y")]), na.rm = TRUE)	

	if(ci.area == TRUE | ci.border == TRUE)
    {
        bounds <- calcResponse(x,alpha=alpha,x.levels=xd)
		bounds.add <- calcResponse(x,alpha=alpha,x.levels=xd.add)
        
		if(equal.axis==TRUE)
        {
			xd <- seq(tmp.range[1],tmp.range[2],length.out=xn)
			xd <- union(xd, tmp.range)
						
			## additional points outer range for nice bounds
			delta <- abs(rx[1]-rx[2])/xn
			xd <- xd[order(xd)]
			xd.add <- c(xd[1]-delta*1:10, xd, xd[length(xd)]+delta*1:10)
			
			bounds <- calcResponse(x,alpha=alpha,x.levels=xd)
			bounds.add <- calcResponse(x,alpha=alpha,x.levels=xd.add)
			
			yrange <- range(c(as.vector(x@data[,c("x","y")]),
                            as.vector(bounds[,c("X","Y","Y.LCI","Y.UCI")])),
                            na.rm=TRUE)		
		}
        else
        {
            yrange <- range(c(as.vector(x@data[,"y"]),
                            as.vector(bounds[,c("Y","Y.LCI","Y.UCI")])),
                            na.rm=TRUE)
        }
    } # end if(ci.area == TRUE | ci.border == TRUE) 
    else
    {
        if(equal.axis==TRUE)
            yrange <- tmp.range
        else
            yrange <- range(as.vector(x@data[,"y"]),na.rm=TRUE)
    }    
                     
    if(equal.axis==TRUE)
    {
        if(is.null(ylim)) xlim <- ylim <- tmp.range
        else xlim <- ylim
  	}
    else
    {
        if(is.null(xlim)) xlim <- rx
        if(is.null(ylim)) ylim <- yrange
    }		
	
	if(is.null(main)) 
        main <- paste(titname,"Fit")
	   
	if(!add)
        plot(0,0,cex=0,ylim=ylim,xlim=xlim,xlab=x.lab,ylab=y.lab,main = main, sub="", bty="n", ...)
        
    else
    {
        sub <- ""
        add.legend <- FALSE
		add.grid <- FALSE
    }
  
 	if(add.cor == TRUE)
    {
        cor.coef <- paste(round(cor(x@data[,"x"],x@data[,"y"],use="pairwise.complete.obs", method=cor.method),3))
        if (cor.method == "pearson") cortext<- paste("Pearson's r = ",cor.coef, sep="")
        if (cor.method == "kendall") cortext<- bquote(paste("Kendall's ", tau, " = ", .(cor.coef) , sep=""))
        if (cor.method == "spearman") cortext <- bquote(paste("Spearman's ", rho, " = ",.(cor.coef), sep=""))
        mtext(side=1,line=-2,cortext,adj=0.9,font=1)
	}
    
    if(ci.area == TRUE | ci.border == TRUE)
    {
        if(ci.area == TRUE)
        {
            xxx <- c(xd.add,xd.add[order(xd.add,decreasing=TRUE)])
            yy1<-c(as.vector(bounds.add[,"Y.LCI"]))
            yy2<-c(as.vector(bounds.add[,"Y.UCI"]))
            yyy <-c(yy1,yy2[order(xd.add,decreasing=TRUE)])
            polygon(xxx,yyy,col=ci.area.col,border="white", lty=0)
        } 
		
		if(add.grid) grid()
		
        if(ci.border == TRUE)
        {
            points(xd.add,bounds.add[,"Y.LCI"], lty=ci.border.lty, lwd=ci.border.lwd, type="l", col=ci.border.col)
            points(xd.add,bounds.add[,"Y.UCI"], lty=ci.border.lty, lwd=ci.border.lwd, type="l", col=ci.border.col)
        }
       
        if(is.null(sub))
        {
            if(x@cimeth %in% c("bootstrap","nestedbootstrap"))
                subtext <- paste("The ", 1-x@alpha,"-confidence bounds are calculated with the ",x@cimeth,"(",x@bootcimeth,") method.",sep="")
            else if((x@regmeth=="PaBa")&(x@cimeth=="analytical"))
                subtext <- ""
            else 
                subtext <- paste("The ", 1-x@alpha,"-confidence bounds are calculated with the ",x@cimeth," method.",sep="")
        }
        else subtext <- sub
    }
    else
    {
		if(add.grid) grid()
        subtext <- ifelse(is.null(sub),"",sub)
    }
   
    if(draw.points == TRUE)  
        points(x@data[,2:3], col=points.col, pch=points.pch, cex=points.cex, ...)
    title(sub=subtext)

#-
  
    if(reg == TRUE)
    {
        b0 <- x@para["Intercept","EST"]
        b1 <- x@para["Slope","EST"]
        abline(b0, b1,,lty=reg.lty,lwd=reg.lwd,col=reg.col )
	}
    
    if(identity == TRUE)   abline(0,1,lty=identity.lty,lwd=identity.lwd, col=identity.col)    
	if(add.legend == TRUE)
    {
        if(identity==TRUE & reg==TRUE)
        {
            text2 <- "identity"
            text1 <- paste(round(x@para["Intercept","EST"],2)," + ",round(x@para["Slope","EST"],2)," * ",x@mnames[1],sep="")
            legend(legend.place,lwd=c(reg.lwd,identity.lwd),lty=c(reg.lty,identity.lty),col=c(reg.col,identity.col),
                   title=paste(titname,"Fit (n=",dim(x@data)[1],")", sep=""),
                   legend=c(text1,text2),box.lty="blank",cex=0.8,bg="white", inset=c(0.01,0.01))
        }
    
        if(identity==TRUE & reg==FALSE)
        {
            text2 <- "identity"
            legend(legend.place,lwd=identity.lwd,lty=identity.lty,col=identity.col,
                    title=paste(titname,"Fit (n=",dim(x@data)[1],")", sep=""),
        	          legend=text2,box.lty="blank",cex=0.8,bg="white", inset=c(0.01,0.01))
        }
    
        if(identity==FALSE & reg==TRUE)
        {
            text1 <- paste(round(x@para["Intercept","EST"],2),"+",
                           round(x@para["Slope","EST"],2),"*",x@mnames[1],sep="")
            legend(legend.place,lwd=c(2),col=c(reg.col),
                   title=paste(titname,"Fit (n=",dim(x@data)[1],")", sep=""),
                   legend=c(text1),box.lty="blank",cex=0.8,bg="white",inset=c(0.01,0.01))
        }
    }
 
	box()
}                                                       


#' Plot Estimated Systematical Bias with Confidence Bounds
#' 
#' This function plots the estimated systematical bias 
#' \eqn{( Intercept + Slope * Refrencemethod ) - Referencemethod}
#' with confidence bounds, covering the whole range of reference method X 
#' or only part of it.
#'
#' @param x object of class "MCResult".
#' @param xn # number of poits for drawing of confidence bounds/area. 
#' @param add logical value. If \code{add=TRUE}, the grafic will be drawn in current grafical window.
#' @param prop a logical value. If \code{prop=TRUE} the proportional bias \eqn{ % bias(Xc) = [ Intercept + (Slope-1) * Xc ] / Xc} will be drawn.
#' @param xlim limits of the x-axis. If \code{xlim=NULL} the x-limits will be calculated automatically.
#' @param ylim limits of the y-axis. If \code{ylim=NULL} the y-limits will be calculated automatically.
#' @param bias logical value. If \code{identity=TRUE} the bias line will be drawn. If ci.bounds=FALSE and ci.area=FALSE the bias line will be drawn always.
#' @param bias.col color of the bias line.
#' @param bias.lty type of the bias line.
#' @param bias.lwd width of the bias line.
#' @param zeroline logical value. If \code{zeroline=TRUE} the zero-line will be drawn.
#' @param zeroline.col color of the zero-line.
#' @param zeroline.lty type of the zero-line.
#' @param zeroline.lwd width of the zero-line.
#' @param ci.area logical value. If \code{ci.area=TRUE} (default) the confidence area will be drawn.
#' @param ci.border logical value. If ci.border=TRUE  the confidence limits will be drawn.
#' @param ci.area.col color of the confidence area.  
#' @param ci.border.col color of the confidence limits. 
#' @param ci.border.lty line type of confidence limits. 
#' @param ci.border.lwd line width of confidence limits.
#' @param cut.point  numeric value. Decision level of interest.
#' @param cut.point.col color of the confidence bounds at the required decision level.
#' @param cut.point.lty line type of the confidence bounds at the required decision level.
#' @param cut.point.lwd line width of the confidence bounds at the required decision level.
#' @param main character string. The main title of plot. If \code{main = NULL} it will include regression name.
#' @param sub character string. The subtitle of plot. If \code{sub=NULL} and \code{ci.border=TRUE} or \code{ci.area=TRUE} it will include the art of confidence bounds calculation.
#' @param add.grid logical value. If \code{grid=TRUE} (default) the gridlines will be drawn.  
#' @param xlab label for the x-axis
#' @param ylab label for the y-axis
#' @param alpha numeric value specifying the 100(1-\code{alpha})\% confidence level of confidence intervals (Default is 0.05).
#' @param ... further graphical parameters
#' @seealso \code{\link{calcBias}}, \code{\link{plot.mcr}}, \code{\link{plotResiduals}}, \code{\link{plotDifference}}, \code{\link{compareFit}}
#' @aliases plotBias
#' @examples
#' #library("mcr")
#' data(creatinine,package="mcr")
#' 
#' creatinine <- creatinine[complete.cases(creatinine),]
#' x <- creatinine$serum.crea
#' y <- creatinine$plasma.crea
#' 
#' # Calculation of models
#' m1 <- mcreg(x,y,method.reg="WDeming", method.ci="jackknife",
#'                 mref.name="serum.crea",mtest.name="plasma.crea", na.rm=TRUE)
#' m2 <- mcreg(x,y,method.reg="WDeming", method.ci="bootstrap",
#'                 method.bootstrap.ci="BCa",mref.name="serum.crea",
#'                 mtest.name="plasma.crea", na.rm=TRUE)
#' 
#' # Grafical comparison of systematical Bias of two models
#' plotBias(m1, zeroline=TRUE,zeroline.col="black",zeroline.lty=1,
#'                 ci.area=TRUE,ci.border=FALSE, ci.area.col=grey(0.9),
#'                 main = "Bias between serum and plasma creatinine",
#'                 sub="Comparison of Jackknife and BCa-Bootstrap confidence bounds ")
#' plotBias(m2, ci.area=FALSE, ci.border=TRUE, ci.border.lwd=2,
#'                 ci.border.col="red",bias=FALSE ,add=TRUE)
#' includeLegend(place="topleft",models=list(m1,m2), lwd=c(10,2),
#'                 lty=c(2,1),colors=c(grey(0.9),"red"), bias=TRUE,
#'                 design="1", digits=4)
#' 
#' # Drawing of proportional bias
#' plotBias(m1, ci.area=FALSE, ci.border=TRUE)
#' plotBias(m1, ci.area=FALSE, ci.border=TRUE, prop=TRUE)
#' plotBias(m1, ci.area=FALSE, ci.border=TRUE, prop=TRUE, cut.point=0.6)
#' plotBias(m1, ci.area=FALSE, ci.border=TRUE, prop=TRUE, cut.point=0.6,
#'              xlim=c(0.4,0.8),cut.point.col="orange", cut.point.lwd=3, main ="")
MCResult.plotBias<-function(x, 
                            xn = 100, 
                            alpha = 0.05, 
                            add = FALSE,
                            prop = FALSE,
                            xlim =NULL, 
                            ylim=NULL,
                            bias = TRUE,
                            bias.lty = 1,
                            bias.lwd = 2,
                            bias.col = NULL,
                            ci.area = TRUE, 
                            ci.area.col = NULL,
                            ci.border = FALSE, 
                            ci.border.col = NULL, 
                            ci.border.lwd = 1,
                            ci.border.lty = 2,
                            zeroline = TRUE,
                            zeroline.col = NULL,
                            zeroline.lty = 2,
                            zeroline.lwd = 1,
                            main = NULL,
                            sub = NULL,
                            add.grid= TRUE,
                            xlab = NULL,
                            ylab = NULL,
                            cut.point = NULL,
                            cut.point.col = "red",
                            cut.point.lwd = 2,
                            cut.point.lty = 1,
                            ...)
{                  
	stopifnot(is.logical(ci.area))
	stopifnot(is.logical(bias))
	stopifnot(is.logical(add))
	stopifnot(is.logical(ci.border))
	stopifnot(is.logical(add.grid))
	stopifnot(is.numeric(alpha))
	stopifnot(is.logical(prop))
	stopifnot(alpha>0 & alpha<1)
	if (!is.null(cut.point)) 
    {
       stopifnot(is.numeric(cut.point))
	    }
    stopifnot(is.numeric(xn))
    stopifnot(round(xn) == xn)
    stopifnot(xn>2)
    if(x@regmeth == "LinReg")
    {
        titname<-"Linear Regression"
	} else if(x@regmeth == "WLinReg")
    {
        titname<-"Weighted Linear Regression"
    } else if(x@regmeth == "Deming")
    {
        titname<-"Deming Regression"
    }  else if(x@regmeth == "WDeming")
    {
        titname<-"Weighted Deming Regression"
    }  else   titname <- "Passing Bablok Regression"
    
	if(x@regmeth %in% c("WDeming","PaBa", "PaBaLarge") & x@cimeth== "analytical" & (ci.area==TRUE | ci.border==TRUE))
    {
        ci.area <- FALSE
        ci.border <- FALSE
        cat("Calculation of confidence bounds with analytical method\n is not available for this regression type.\n")
    }
 
    if(is.null(xlab))  xlab <- x@mnames[1]
  	if(is.null(ylab))  if ( prop == FALSE ) ylab <- "Bias" else ylab <- "% Bias"
  	if (!is.null(xlim) & !is.null(cut.point)) 
    {
          cut.point <- cut.point[cut.point<=xlim[2] & cut.point>=xlim[1] ]
          if (length(cut.point)==0) cut.point <- NULL
	}
    
    b0<-x@para["Intercept","EST"]
    b1<-x@para["Slope","EST"]
    b0.low <- x@para["Intercept","LCI"]
    
    stopifnot( b0.low < b0 )
   
    niceblue <- rgb(37/255,52/255,148/255)
    niceor <- rgb(230/255,85/255,13/255)
    niceblue.bounds <- rgb(236/255, 231/255, 242/255)
    
    if(is.null(bias.col)) 
          bias.col <- niceblue  
  	if(is.null(ci.area.col)) 
          ci.area.col <- niceblue.bounds
  	if(is.null(ci.border.col)) 
          ci.border.col <- niceblue
  	if(is.null(zeroline.col)) 
          zeroline.col <- "black"
          
    if(is.null(sub))
    {
        if(x@cimeth %in% c("bootstrap","nestedbootstrap"))
            subtext <- paste("The ", 1-x@alpha,"-confidence bounds are calculated with the ",x@cimeth,"(",x@bootcimeth,") method.",sep="")
        else if((x@regmeth %in% c("PaBa", "PaBaLarge"))&(x@cimeth=="analytical"))
            subtext <- ""
        else 
            subtext <- paste("The ", 1-x@alpha,"-confidence bounds are calculated with the ",x@cimeth," method.",sep="") 
    } else subtext <- sub
           
  #---------- ordinar bias -----------------------------------------------------
          
    if (prop == FALSE)
    {
        if(is.null(main)) main <- paste("Bias plot \n",titname,"Fit")
        if (is.null(xlim))
        {
            xlim <-  range(x@data[,"x"],na.rm=TRUE)
            if (!is.null(cut.point))
            {
                if ( min(cut.point) <=xlim[1] & min(cut.point) < 0 ) xlim[1] <- 1.1* min(cut.point)
                if ( min(cut.point) <=xlim[1] &  min(cut.point) > 0 ) xlim[1] <- 0.9* min(cut.point)
                if ( max(cut.point) >=xlim[2] & max(cut.point) < 0 ) xlim[2] <- 0.9*max(cut.point)
                if ( max(cut.point) >=xlim[2] & max(cut.point)> 0 ) xlim[2] <- 1.1*max(cut.point)
                if ( max(cut.point)>=xlim[2] & max(cut.point) == 0 ) xlim[2] <- max(cut.point)+0.1*abs(xlim[2]-xlim[1])
                if ( min(cut.point)<=xlim[1] &  min(cut.point) == 0 ) xlim[1] <-  min(cut.point)-0.1*abs(xlim[2]-xlim[1])
            }    
        }
        
        ## points for drawing of confidence bounds/area:
        xd <- seq(xlim[1],xlim[2],length.out=xn)
        #if (xlim[1]==0) xd <- seq(0,xlim[2]+0.05*(xlim[2]-xlim[1]),length.out=xn)
        xd <- union(union(xd, xlim), cut.point)
        xd <- xd[order(xd)]
       
    
        bounds<-calcBias(x,x.levels=xd,alpha=alpha, type="absolute")
    
        if(ci.area == TRUE | ci.border == TRUE)
        {	
            if(is.null(ylim)==FALSE & add==FALSE)
            {
                yrange <- ylim
            } else  yrange<-range(bounds[,c("Bias","LCI","UCI")],na.rm=TRUE)
	
    	    if(add == FALSE) 
            {
                plot(0,0,cex=0, xlim=range(xd), ylim=yrange, xlab=xlab, 
                     ylab=ylab, main = main, sub="", bty="n", ...)
            } 
            else
            {
                subtext <- ""
                legend <- FALSE
                grid <- FALSE
                zeroline <- FALSE
            }
            if(ci.area == TRUE)
            {
                xxx <- c(xd,xd[order(xd,decreasing=TRUE)])
                yy1<-c(as.vector(bounds[,"LCI"]))
                yy2<-c(as.vector(bounds[,"UCI"]))
                yyy <-c(yy1,yy2[order(xd,decreasing=TRUE)])
                polygon(xxx,yyy,col=ci.area.col,border="white")
                if(add.grid== TRUE) grid() 
            } #end if(ci.area == TRUE)
            
            if(ci.border == TRUE)
            {
                points(xd, bounds[,"LCI"], lty=ci.border.lty,lwd=ci.border.lwd, col=ci.border.col, type="l")
          	    points(xd, bounds[,"UCI"], lty=ci.border.lty,lwd=ci.border.lwd, col=ci.border.col, type="l")
            } #end if(ci.border == TRUE)
            if(bias == TRUE) points(xd, bounds[,"Bias"], col=bias.col, type="l", lwd=bias.lwd, lty=bias.lty)
            if(zeroline == TRUE) abline(h=0,lty=zeroline.lty, col=zeroline.col, lwd=zeroline.lwd)
            title(sub = subtext)  
        } else {
            if(zeroline == TRUE & is.null(ylim) & is.null(cut.point))   ylim <- range(c(b0+b1*xd-xd, 0)) 
            if(zeroline == TRUE & is.null(ylim) & !is.null(cut.point))  ylim <- range(c(b0+b1*xd-xd, 0,
                                                                           calcBias(x, x.levels=cut.point, alpha=alpha)[,c("Bias","LCI","UCI")])) 
            if(zeroline == FALSE & is.null(ylim) & is.null(cut.point))  ylim <- range(c(b0+b1*xd-xd)) 
            if(zeroline == FALSE & is.null(ylim) & !is.null(cut.point)) ylim <- range(c(b0+b1*xd-xd, 
                                                                           calcBias(x, x.levels=cut.point, alpha=alpha)[,c("Bias","LCI","UCI")])) 
      
            if (add == FALSE)
            {
                plot(xd,b0+b1*xd-xd, main=main, col=bias.col, lty=bias.lty, ylim=ylim, 
                     lwd=bias.lwd, type="l",xlab=xlab,ylab=ylab,sub="",...)
            } else {
                main <- ""
                subtext <- ""
                add.grid<- FALSE
                zeroline <- FALSE
                points(xd,b0+b1*xd-xd, main=main, col=bias.col, lty=bias.lty, ylim=ylim, 
                       lwd=bias.lwd, type="l",xlab=xlab,ylab=ylab,sub="",...)
            } #end if (add=FALSE)
        }#end if(ci.area == TRUE | ci.border == TRUE)
    
        title(sub = subtext)       
        if(add.grid== TRUE) grid() 
        if(zeroline == TRUE) abline(h=0,lty=zeroline.lty, col=zeroline.col, lwd=zeroline.lwd)    
        if (!is.null(cut.point))
        {
            cut.point.y <- calcBias(x,x.levels=cut.point,alpha=alpha)
            for (ii in seq_along(cut.point)){
                     segments(cut.point[ii], cut.point.y[ii,"LCI"],cut.point[ii], cut.point.y[ii,"UCI"],
                     col=cut.point.col, lty=cut.point.lty,lwd=cut.point.lwd)
                     }
        }
    }
    
    # --------  proportional bias ------------------------------------------------
    if (prop == TRUE)
    {   
        ## This function calculates x-points for drawing of confidence bounds/area
        golden.seq <- function(x.low, x.up,xn, result)
        {
            delta <- abs(x.up-x.low)/1000000
            p <- abs(1-(delta/(x.up-x.low))^(1/xn))
            
            temp <- function(x.low, x.up,p,xn,tresh, result)
            {
                if (abs(x.up-x.low)>tresh)
                {
                    x.up <-(1-p)*x.up+p*x.low
                    result <-temp(x.low, x.up,p,xn,tresh, result)
                    return(c(result,x.up))
                }
            }
            return(temp(x.low, x.up,p,xn,tresh=delta, result))
        } # end golden.seq
    
     
        if(is.null(main))   main <- paste("Proportional bias plot \n",titname,"Fit")
        if (!is.null(cut.point) )
        {
            if (any(cut.point == 0))
            { 
                warning("The proportional bias  at 0 is infinite.")
                cut.point <- NULL
            }
        } 
     
        c <- 0.4   # This coefficient controls the default x-limits of the graph.
      
        if (b0==0)
        {
            lowlimit <- c*sqrt(abs(b0.low))
        } else {
            lowlimit <- c*sqrt(abs(b0))
        }  
     
        if (is.null(xlim))
        {
            xlim <- range(x@data[,"x"])
            if (!is.null(cut.point))
            {
                if ( min(cut.point) <=xlim[1] & min(cut.point) < 0 ) xlim[1] <- 1.1* min(cut.point)
                if ( min(cut.point) <=xlim[1] &  min(cut.point) > 0 ) xlim[1] <- 0.9* min(cut.point)
                if ( max(cut.point) >=xlim[2] & max(cut.point) < 0 ) xlim[2] <- 0.9*max(cut.point)
                if ( max(cut.point) >=xlim[2] & max(cut.point)> 0 ) xlim[2] <- 1.1*max(cut.point)
            }    
        }   
   
        xlim.new <- xlim
    
        if ( xlim[1]==0 )  xlim.new[1] <- min(lowlimit,xlim[2],na.rm=TRUE)/2
        if ( xlim[2]==0 )  xlim.new[2] <- -1*min(lowlimit,abs(xlim[1]),na.rm=TRUE)/2
    
        ## in this case the positive and negative part of the graph 
        ## must be drawn separately
        if (xlim.new[1]<0 & xlim.new[2]>0)
        {
            lowlimit <-min(abs(xlim.new[1])/5,abs(xlim.new[2])/5,lowlimit/5)       
            xlim.neg <- c(xlim.new[1], -lowlimit)
            xd.neg <- -1*golden.seq(x.low=min(abs(xlim.neg[1]),abs(xlim.neg[2])), x.up=max(abs(xlim.neg[1]),abs(xlim.neg[2])),
                                    xn=xn-1,result=max(abs(xlim.neg[1]),abs(xlim.neg[2]))) 
            xd.neg <- union(union(xd.neg, xlim.neg), cut.point[cut.point<0])
            xd.neg <-xd.neg[order(xd.neg)]
            xd.neg<- xd.neg[xd.neg!=0]
            bounds.neg<-calcBias(x,x.levels=xd.neg,alpha=alpha, type="proportional", percent=TRUE) 
          
            xlim.pos <- c( lowlimit , xlim.new[2] )
            xd.pos <- golden.seq(x.low=min(abs(xlim.pos[1]),abs(xlim.pos[2])), x.up=max(abs(xlim.pos[1]),abs(xlim.pos[2])),
                                 xn=xn-1, result=max(abs(xlim.pos[1]),abs(xlim.pos[2])))  
            xd.pos <- union(union(xd.pos, xlim.pos), cut.point[cut.point>0])
            xd.pos<- xd.pos[order(xd.pos)]
            xd.pos<- xd.pos[xd.pos!=0]
            bounds.pos<-calcBias(x,x.levels=xd.pos,alpha=alpha, type="proportional", percent=TRUE) 
            xd <- c(xd.neg, xd.pos)
       
        }
        else
        {
            xd <- golden.seq(x.low=min(abs(xlim.new[1]),abs(xlim.new[2])), x.up=max(abs(xlim.new[1]),abs(xlim.new[2])),
                             xn=xn-1, result=max(abs(xlim.new[1]),abs(xlim.new[2]))) 
            xd <- xd*sign(xlim.new[2])
            xd <-  union(union(xd, xlim.new), cut.point)
                 
            xd<- xd[xd!=0]
            xd<- xd[order(xd)]
        }
    
        bounds <- calcBias(x,x.levels=xd,alpha=alpha, type="proportional", percent=TRUE)  
     
        if(ci.area == TRUE | ci.border == TRUE)
        {	  
            if(is.null(ylim)==FALSE & add==FALSE)
            {
                yrange <- ylim
            } else  yrange<-range(bounds[,c("Prop.bias(%)","LCI","UCI")],na.rm=TRUE)
	
    	    if(add == FALSE) 
            {
                plot(0,0,cex=0, xlim=range(xd), ylim=yrange, xlab=xlab, 
                     ylab=ylab, main = main, sub="", bty="n", ...)
            } 
            else
            {
                subtext <- ""
                legend <- FALSE
                grid <- FALSE
                zeroline <- FALSE
            }
            
        #--
    
            if(ci.area == TRUE)
            {
                if (xlim.new[1]<0 & xlim.new[2]>0)
                {  
                    xxx.pos <- c(xd.pos,xd.pos[order(xd.pos,decreasing=TRUE)])
                    yy1.pos<-c(as.vector(bounds.pos[,"LCI"]))
                    yy2.pos<-c(as.vector(bounds.pos[,"UCI"]))
                    yyy.pos <-c(yy1.pos,yy2.pos[order(xd.pos,decreasing=TRUE)])
                    polygon(xxx.pos,yyy.pos,col=ci.area.col,border="white")

                    xxx.neg <- c(xd.neg,xd.neg[order(xd.neg,decreasing=TRUE)])
                    yy1.neg<-c(as.vector(bounds.neg[,"LCI"]))
                    yy2.neg<-c(as.vector(bounds.neg[,"UCI"]))
                    yyy.neg <-c(yy1.neg,yy2.neg[order(xd.neg,decreasing=TRUE)])
                    polygon(xxx.neg,yyy.neg,col=ci.area.col,border="white")      
                }
                else
                {                    
                    xxx <- c(xd,xd[order(xd,decreasing=TRUE)])
                    yy1<-c(as.vector(bounds[,"LCI"]))
                    yy2<-c(as.vector(bounds[,"UCI"]))
                    yyy <-c(yy1,yy2[order(xd,decreasing=TRUE)])
                    polygon(xxx,yyy,col=ci.area.col,border="white")
                } # end if (range(xd,rm.na)[1]<0 & range(xd,rm.na)[2]>0)
    
                if(add.grid== TRUE)  grid() 
            } #end if(ci.area == TRUE)
        
            if(ci.border == TRUE)
            {
                if (xlim.new[1]<0 & xlim.new[2]>0)
                {
                    points(xd.pos, bounds.pos[,"LCI"], lty=ci.border.lty,lwd=ci.border.lwd, col=ci.border.col, type="l")
      	            points(xd.neg, bounds.neg[,"UCI"], lty=ci.border.lty,lwd=ci.border.lwd, col=ci.border.col, type="l")
                    points(xd.pos, bounds.pos[,"UCI"], lty=ci.border.lty,lwd=ci.border.lwd, col=ci.border.col, type="l")
      	            points(xd.neg, bounds.neg[,"LCI"], lty=ci.border.lty,lwd=ci.border.lwd, col=ci.border.col, type="l")
                }
                else
                {
                    points(xd, bounds[,"LCI"], lty=ci.border.lty,lwd=ci.border.lwd, col=ci.border.col, type="l")
                    points(xd, bounds[,"UCI"], lty=ci.border.lty,lwd=ci.border.lwd, col=ci.border.col, type="l")
                } #end if (range(xd,rm.na)[1]<0 & range(xd,rm.na)[2]>0)
            } #end if(ci.border == TRUE)
            
            if(bias == TRUE) 
            {
                if (xlim.new[1]<0 & xlim.new[2]>0)
                {
                    points(xd.pos, bounds.pos[,"Prop.bias(%)"], col=bias.col, type="l", lwd=bias.lwd, lty=bias.lty)
                    points(xd.neg, bounds.neg[,"Prop.bias(%)"], col=bias.col, type="l", lwd=bias.lwd, lty=bias.lty)
                }
                else  points(xd, bounds[,"Prop.bias(%)"], col=bias.col, type="l", lwd=bias.lwd, lty=bias.lty)   
            } # end  if(bias == TRUE)
            
        }
        else
        { 
            if(zeroline == TRUE & is.null(ylim) & is.null(cut.point))   ylim <- range(c(100*(b0+b1*xd-xd)/xd, 0)) 
            if(zeroline == TRUE & is.null(ylim) & !is.null(cut.point))  ylim <- range(c(100*(b0+b1*xd-xd)/xd, 0,
                                                                           calcBias(x, x.levels=cut.point, 
                                                                           alpha=alpha, type="proportional", 
                                                                           percent=TRUE)[,c("Prop.bias(%)","LCI","UCI")])) 
            if(zeroline == FALSE & is.null(ylim) & is.null(cut.point))  ylim <- range(100*(b0+b1*xd-xd)/xd) 
            if(zeroline == FALSE & is.null(ylim) & !is.null(cut.point)) ylim <- range(c(100*(b0+b1*xd-xd)/xd,
                                                                           calcBias(x, x.levels=cut.point, 
                                                                           alpha=alpha, type="proportional", 
                                                                           percent=TRUE)[,c("Prop.bias(%)","LCI","UCI")]))
            if (add == FALSE)
            {
                if (xlim.new[1]<0 & xlim.new[2]>0)
                {
                    plot(xd.pos,100*(b0+b1*xd.pos-xd.pos)/xd.pos, main=main, col=bias.col, lty=bias.lty, ylim=ylim, 
                         lwd=bias.lwd, type="l",xlab=xlab,ylab=ylab,sub="",xlim=xlim.new,...)
                    points(xd.neg,100*(b0+b1*xd.neg-xd.neg)/xd.neg,col=bias.col, lty=bias.lty,  
                           lwd=bias.lwd, type="l")
                }
                else  plot(xd,100*(b0+b1*xd-xd)/xd, main=main, col=bias.col, lty=bias.lty, ylim=ylim, 
                           lwd=bias.lwd, type="l",xlab=xlab,ylab=ylab,sub="",...)
            }
            else
            {
                main <- ""
                subtext <- ""
                add.grid<- FALSE
                zeroline <- FALSE
                if (xlim[1]<0 & xlim[2]>0)
                {
                    points(xd.pos,100*(b0+b1*xd.pos-xd.pos)/xd.pos, main=main, col=bias.col, 
                    lty=bias.lty, ylim=ylim,lwd=bias.lwd, type="l",xlab=xlab,ylab=ylab,sub="",...)
                    points(xd.neg,100*(b0+b1*xd.neg-xd.neg)/xd.neg, main=main, col=bias.col, 
                    lty=bias.lty, ylim=ylim, lwd=bias.lwd, type="l",xlab=xlab,ylab=ylab,sub="",...)
                }
                else points(xd,100*(b0+b1*xd-xd)/xd, main=main, col=bias.col, 
                     lty=bias.lty, ylim=ylim,lwd=bias.lwd, type="l",xlab=xlab,ylab=ylab,sub="",...)
            }#end if (add=FALSE)
        } #end if(ci.area == TRUE | ci.border == TRUE)
                                 
        title(sub = subtext)       
        if(add.grid== TRUE) grid() 
        if(zeroline == TRUE) abline(h=0,lty=zeroline.lty, col=zeroline.col, lwd=zeroline.lwd)   
        if (!is.null(cut.point))
        {
            cut.point.y <- calcBias(x,x.levels=cut.point,alpha=alpha, type="proportional", percent=TRUE)
            for (ii in seq_along(cut.point)){
                     segments(cut.point[ii], cut.point.y[ii,"LCI"],cut.point[ii], cut.point.y[ii,"UCI"],
                     col=cut.point.col, lty=cut.point.lty,lwd=cut.point.lwd)
                     }
        }          
        if (xlim[1]<0 & xlim[2]>0) abline(v=0, lty=zeroline.lty, col=zeroline.col)
    } # end if (prop == TRUE)
    
    box()
}       


#' Print Summary of a Regression Analysis
#'
#' @param .Object object of type "MCResult".
#' @seealso \code{\link{getCoefficients}}, \code{\link{getRegmethod}}
#' @aliases printSummary
MCResult.printSummary<-function(.Object)
{
    print("MCResult virtual function")
    return(NULL)
}

#' Plot Residuals of an MCResult Object
#'
#' @param .Object object of type "MCResult".
#' @param res.type If res.type="y" the difference between the test method and it's prediction will be drawn. 
#' If res.type="x" the reference method and it's prediction will be drawn. 
#' In case ordinary and weighted ordinary linear regression 
#' this difference will be zero.  
#' @param xaxis Values on the x-axis. One can choose from estimated values of x (xaxis=\code{"xhat"}), 
#' y (xaxis=\code{"xhat"}) or the mean of estimated values of x and y (\code{xaxis="both"}).
#' If res.type="optimized" the proper type of residuals for each regression will be drawn. 
#' @param ref.line logical value. If \code{ref.line = TRUE} (default), the reference line will be drawn.
#' @param ref.line.col reference line color.
#' @param ref.line.lty reference line type.
#' @param ref.line.lwd reference line width.
#' @param xlab label for the x-axis
#' @param ylab label for the y-axis
#' @param add.grid logical value. If \code{add.grid = TRUE} (default) the gridlines will be drawn.
#' @param main character string specifying the main title of the plot
#' @param ... further graphical parameters
#' @seealso \code{\link{getResiduals}}, \code{\link{plot.mcr}}, \code{\link{plotDifference}}, \code{\link{plotBias}}, \code{\link{compareFit}}
#' @aliases plotResiduals
#' @examples
#'     data(creatinine,package="mcr")
#'     x <- creatinine$serum.crea
#'     y <- creatinine$plasma.crea
#' 
#'     # Deming regression fit.
#'     # The confidence intercals for regression coefficients
#'     # are calculated with analytical method
#'     model <- mcreg( x,y,error.ratio=1,method.reg="WDeming", method.ci="jackknife",
#'                      mref.name = "serum.crea", mtest.name = "plasma.crea", na.rm=TRUE )
#'     plotResiduals(model, res.type="optimized", xaxis="both" )
#'     plotResiduals(model, res.type="y", xaxis="yhat")
MCResult.plotResiduals<-function(.Object, res.type=c("optimized", "y", "x"),
                                    xaxis=c("yhat","both","xhat"),
                                    ref.line=TRUE,
                                    ref.line.col="red",
                                    ref.line.lty=2,
                                    ref.line.lwd=1,
                                    main=NULL,
                                    xlab=NULL,
                                    ylab=NULL,
                                    add.grid=TRUE,
                                    ...)
{    
    res.type <- match.arg(res.type)
    xaxis <- match.arg(xaxis)
    res <- getResiduals(.Object)
    r <- res[,res.type]
    if (xaxis=="yhat") est <- getData(.Object)[,"Y"]- res[,"y"]
    if (xaxis=="xhat") est <- getData(.Object)[,"X"]- res[,"x"]
    if (xaxis=="both") est <- ((getData(.Object)[,"X"]- res[,"x"])+(getData(.Object)[,"Y"]- res[,"y"]))/2
    
    if(is.null(main))
    {
        if(.Object@regmeth == "LinReg")
		    titname <- "Linear Regression"
        else if(.Object@regmeth == "WLinReg")
		    titname <- "Weighted Linear Regression"
        else if(.Object@regmeth == "Deming")
		    titname <- "Deming Regression"
        else if(.Object@regmeth == "WDeming")
            titname <- "Weighted Deming Regression"
        else
            titname <- "Passing Bablok Regression"

        main <- paste("Residual plot for", titname,"Fit"  ) 
    }
    if(is.null(ylab)) 
    {
        if ( res.type=="x") ylab <- paste("Residuals (",.Object@mnames[1],")", sep="")
        if ( res.type=="y") ylab <- paste("Residuals (",.Object@mnames[2],")", sep="")
        if ( res.type=="optimized") ylab <- "Optimized residuals"
    } else{}    
    if(is.null(xlab))
    {     
        if (xaxis=="yhat") xlab <- paste("Estimated values of",.Object@mnames[2])
        if (xaxis=="xhat") xlab <- paste("Estimated values of",.Object@mnames[1])
        if (xaxis=="both") xlab <- paste("Mean of estimated values of both methods")   
    }
    plot(est,r,main=main,xlab=xlab,ylab=ylab, ...)
    
    if(add.grid == TRUE) 
        grid()
    if(ref.line == TRUE) 
        abline(h=0, col=ref.line.col, lty=ref.line.lty, lwd=ref.line.lwd)
}

