###############################################################################
##
## MCResultBCaMethods.R
##
## Definition of methods for class MCResultBCa
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

#' MCResultBCa object constructor with matrix in wide format as input.
#'
#' @param wdata Measurement data in matrix format. First column reference method (x), second column comparator method (y).
#' @param para Regression parameters in matrix form. Rows: Intercept, Slope. Cols: EST, SE, LCI, UCI.
#' @param xmean Global (weighted) mean of x-values
#' @param sample.names Names of individual data points, e.g. barcodes of measured samples.
#' @param method.names Names of reference and comparator method.
#' @param regmeth Name of statistical method used for regression.
#' @param cimeth Name of statistical method used for computing confidence intervals.
#' @param error.ratio Ratio between standard deviation of reference and comparator method.
#' @param alpha 1 - significance level for confidence intervals.
#' @param glob.coef Numeric vector of length two with global point estimations of intercept and slope.
#' @param glob.sigma Numeric vector of length two with global estimations of standard errors of intercept and slope.
#' @param bootcimeth Bootstrap based confidence interval estimation method.
#' @param nsamples Number of bootstrap samples.
#' @param nnested Number of nested bootstrap samples.
#' @param rng.seed Seed used to call mcreg, NULL if no seed was used
#' @param rng.kind RNG type (string, see set.seed for details) used, only meaningfull if rng.seed was specified
#' @param B0jack Numeric vector with point estimations of intercept for jackknife samples.
#' @param B1jack Numeric vector with point estimations of slope for jackknife samples.
#' @param B0 Numeric vector with point estimations of intercept for each bootstrap sample.
#' @param B1 Numeric vector with point estimations of slope for each bootstrap sample.
#' @param sigmaB0 Numeric vector with estimation of standard error of intercept for each bootstrap sample.
#' @param sigmaB1 Numeric vector with estimation of standard error of slope for each bootstrap sample.
#' @param MX Numeric vector with point estimations of (weighted-)average of reference method values for each bootstrap sample. 
#' @param weight numeric vector specifying the weights used for each point
#' @return MCResult object containing regression results.
newMCResultBCa <- function( wdata, para, xmean, sample.names=NULL, method.names=NULL, regmeth="unknown", glob.coef, glob.sigma,
                            cimeth="unknown", bootcimeth="unknown", nsamples, nnested, rng.seed,rng.kind, B0jack, B1jack,
                            B0, B1, MX, sigmaB0, sigmaB1, error.ratio, alpha=0.05,
                            weight=rep(1,nrow(wdata)))  
{
    ## Check validity of parameters
    stopifnot(is.numeric(xmean))
    stopifnot(is.matrix(wdata))
    stopifnot(ncol(wdata)==2)
    stopifnot(is.matrix(para))
    stopifnot(all(dim(para)==c(2,4)))
    stopifnot(is.character(regmeth))
    stopifnot(is.element(regmeth,c("LinReg","WLinReg","Deming","PaBa","WDeming", "PaBaLarge")))
    stopifnot(is.character(cimeth))
    stopifnot(is.element(cimeth,c("bootstrap","nestedbootstrap")))
    stopifnot(is.character(bootcimeth))
    stopifnot(is.element(bootcimeth,c("BCa")))
    stopifnot(is.numeric(alpha))
    stopifnot(alpha<=1 & alpha>=0)
    stopifnot(is.numeric(glob.coef))
    stopifnot(length(glob.coef)==2)
    stopifnot(is.numeric(glob.sigma))
    stopifnot(length(glob.sigma)==2)
    stopifnot(is.numeric(nsamples))
    stopifnot(round(nsamples)==nsamples)
    stopifnot(nsamples>0)
    stopifnot(is.numeric(nnested))
    stopifnot(round(nnested)==nnested)
    stopifnot(nnested>0)
    
    # RNG settings
    stopifnot(is.numeric(rng.seed))
    stopifnot(is.character(rng.kind))
    
    stopifnot(is.numeric(B0jack))
    stopifnot(is.numeric(B1jack))
    stopifnot(length(B0jack)==length(B1jack))

    stopifnot(is.numeric(B0))
    stopifnot(is.numeric(B1))
    stopifnot(length(B0)==length(B1))
    stopifnot(length(B0)==nsamples)

    stopifnot(is.numeric(MX))
    stopifnot(length(MX)==nsamples)

    stopifnot(is.numeric(sigmaB0))
    stopifnot(is.numeric(sigmaB1))
    stopifnot(length(sigmaB0)==length(sigmaB1))

    stopifnot(error.ratio>=0)
    stopifnot(is.numeric(error.ratio))

    ## Sample names
    if(is.null(sample.names)) 
        snames <- paste("S", 1:nrow(wdata),sep="")
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
    rownames(para) <- c("Intercept", "Slope")
    colnames(para) <- c("EST", "SE", "LCI", "UCI")
    names(weight) <- snames

    new(Class="MCResultBCa", data=data, para=para, xmean=xmean, mnames=mnames, regmeth, cimeth=cimeth, bootcimeth=bootcimeth, alpha=alpha,
        glob.coef=glob.coef, glob.sigma=glob.sigma, nsamples=nsamples, nnested=nnested, rng.seed=rng.seed, rng.kind=rng.kind,
        B0jack=B0jack, B1jack=B1jack, B0=B0, B1=B1, MX=MX, sigmaB0=sigmaB0, sigmaB1=sigmaB1, error.ratio=error.ratio, weight=weight)
}


###############################################################################
## Methods
###############################################################################

#' Initialize Method for 'MCResultBCa' Objects.
#' 
#' Method initializes newly created objects of class 'MCResultBCa'.
#' 
#' @param .Object object to be initialized
#' @param data empty data.frame
#' @param xmean 0 for init-purpose
#' @param para empty coefficient matrix
#' @param mnames empty method names vector
#' @param regmeth string specifying the regression-method 
#' @param cimeth string specifying the confidence interval method
#' @param bootcimeth string specifying the method for bootstrap confidence intervals
#' @param error.ratio for deming regression 
#' @param alpha value specifying the 100(1-alpha)% confidence-level
#' @param glob.coef global coefficients
#' @param rng.seed random number generator seed
#' @param rng.kind type of the random number generator
#' @param glob.sigma global sd values for regression parameters
#' @param nsamples number of samples for resampling
#' @param nnested number of inner simulation for nested bootstrap
#' @param B0jack jackknife intercpet
#' @param B1jack jackknife slope
#' @param B0 intercept
#' @param B1 slope
#' @param MX parameter
#' @param sigmaB0 SD for intercepts
#' @param sigmaB1 SD for slopes
#' @param weight 1 for each data point

MCResultBCa.initialize <- function( .Object, data=data.frame(X=NA,Y=NA), para=matrix(NA,ncol=4,nrow=2), xmean=0,
                                    mnames=c("unknown","unknown"), regmeth="unknown", cimeth="unknown", bootcimeth="unknown", alpha=0.05, glob.coef=c(0,0),
                                    glob.sigma=c(0,0), nsamples=0, nnested=0, B0jack=0, B1jack=0, B0=0, B1=0, MX=0, rng.seed=as.numeric(NA), rng.kind="unknown",
                                    sigmaB0=0, sigmaB1=0, error.ratio=0, weight=1) 
{
    .Object@data <- data
    .Object@error.ratio <- error.ratio
    .Object@para <- para
    .Object@xmean <- xmean
    .Object@mnames <- mnames
    .Object@regmeth <- regmeth
    .Object@alpha <- alpha
    .Object@cimeth<-cimeth
    .Object@bootcimeth<-bootcimeth
    .Object@glob.coef <- glob.coef
    .Object@glob.sigma <- glob.sigma
    .Object@nsamples <- nsamples
    .Object@nnested <- nnested
    .Object@B0jack <- B0jack
    .Object@B1jack <- B1jack
    .Object@B0 <- B0
    .Object@B1 <- B1
    .Object@MX <- MX
    .Object@sigmaB0 <- sigmaB0
    .Object@sigmaB1 <- sigmaB1
    .Object@rng.seed <- rng.seed
    .Object@rng.kind <- rng.kind    
    .Object@weight<-weight    
    return(.Object)
}

#' Plot distriblution of bootstrap coefficients
#'
#' Plot distriblution of bootstrap coefficients (slope and intercept).
#'
#' @param .Object Object of class "MCResultBCa"
#' @param breaks used in function 'hist' (see ?hist)
#' @param ... further graphical parameters

MCResultBCa.plotBootstrapCoefficients<-function(.Object, breaks=20, ...)
{
    layout(matrix(c(1,1,2,2,1,1,3,3),nrow=2,ncol=4, byrow=TRUE))
    plot(.Object@B1,.Object@B0, xlab="bootstrap slope",ylab="bootstrap intercept",...)
    abline(v=.Object@glob.coef[2], col="red",lty=2)
    abline(h=.Object@glob.coef[1], col="red",lty=2)
    dns<-density(.Object@B1)
    hi<-hist(.Object@B1,breaks=breaks,plot=FALSE)
    hist(.Object@B1,freq=FALSE,ylim=range(c(dns$y,hi$density),na.rm=TRUE),breaks=breaks,
            main="Histogram for Slope",xlab="bootstrapped coefficients",border="white",col="lightgrey",ylab="",yaxt="n",...)
    points(density(.Object@B1),type="l",lwd=1.5,lty=2)
    abline(v=.Object@glob.coef[2], col="red")
    text(.Object@glob.coef[2]+(hi$breaks[length(hi$breaks)]-hi$breaks[1])/30,range(c(dns$y,hi$density))[2],"estimation",col="red",adj=0)
    dns<-density(.Object@B0)
    hi<-hist(.Object@B0,breaks=breaks,plot=FALSE)
    hist(.Object@B0,freq=FALSE,ylim=range(c(dns$y,hi$density)),breaks=breaks,
            main="Histogram for Intercept",xlab="bootstrapped coefficients",border="white",ylab="",col="lightgrey",yaxt="n",...)
    points(density(.Object@B0),type="l",lty=2)
    abline(v=.Object@glob.coef[1], col="red")
    text(.Object@glob.coef[1]+(hi$breaks[length(hi$breaks)]-hi$breaks[1])/30,range(c(dns$y,hi$density),na.rm=TRUE)[2],"estimation",col="red",adj=0)
}

#' Plot distriblution of bootstrap pivot T
#'
#' Plot distriblution of bootstrap pivot T for slope and intercept and compare
#' them with t(n-2) distribution.
#'
#' @param .Object Object of class "MCResultBCa".
#' @param breaks Number of breaks in histogram.
#' @param ... further graphical parameters

MCResultBCa.plotBootstrapT<-function(.Object,breaks=20,...)
{
    if (length(.Object@sigmaB0)<=1)
    {
        return("T*-density is not available (there is no analytical SE-estimation for this regression type)")
    } 
    else 
    {
        DF <- length(.Object@data[,"x"])-2
        tstarB0 <- (.Object@B0-.Object@glob.coef[1])/.Object@sigmaB0
        tstarB1 <- (.Object@B1-.Object@glob.coef[2])/.Object@sigmaB1
        RB0 <- range(tstarB0,na.rm=TRUE)
        RB1 <- range(tstarB1,na.rm=TRUE)
        seqRB0 <- seq(RB0[1], RB0[2], by=(RB0[2]-RB0[1])/100)
        seqRB1 <- seq(RB1[1], RB1[2], by=(RB1[2]-RB1[1])/100)
        RT <- rt(100000,df=DF)
        DF <- length(.Object@data[,"x"])-2
        par(mfrow=c(2,3))
        plot(.Object@B1,.Object@sigmaB1, xlab="bootstrap slope", ylab="bootstrap se for slope")
        dns <- density(tstarB1)
        hi <- hist(tstarB1, breaks=breaks, plot=FALSE)
        hist(   tstarB1, freq=FALSE, ylim=range(c(dns$y, hi$density, dt(seqRB1,df=DF)), na.rm=TRUE), breaks=breaks,
                main="t* for slope",xlab="t*",border="white",col="lightgrey",...)
        points(density(tstarB1), type="l", lty=2)
        points(seqRB1, dt(seqRB1, df=DF), type="l", col="red")
        legend("topright", col=c("black", "red"), lty=c(2,1), legend=c("t*", "t(n-2)"), box.lty="blank")
        qqplot( RT, tstarB1, ylab="t*", xlab="Quantiles of t(n-2)", main="t* for slope",
                xlim=range(c(tstarB1, RT), na.rm=TRUE), ylim=range(c(tstarB1,RT), na.rm=TRUE))
        abline(0, 1)
        plot(.Object@B0, .Object@sigmaB0, xlab="bootstrap intercept", ylab="bootstrap se for intercept")
        dns <- density(tstarB0)
        hi <- hist(tstarB0, breaks=breaks, plot=FALSE)
        hist(tstarB0, freq=FALSE, ylim=range(c(dns$y, hi$density, dt(seqRB0, df=DF)), na.rm=TRUE), breaks=breaks,
             main="t* for intercept", xlab="t*", border="white", col="lightgrey", ...)
        points(density(tstarB0), type="l", lty=2)
        points(seqRB0, dt(seqRB0, df=DF), type="l", col="red")
        legend("topright", col=c("black", "red"), lty=c(2,1), legend=c("t*","t(n-2)"), box.lty="blank")
        qqplot( RT, tstarB0, ylab="t*", xlab="Quantiles of t(n-2)", main="t* for intercept", xlim=range(c(tstarB0, RT), na.rm=TRUE),
                ylim=range(c(tstarB0, RT), na.rm=TRUE))
        abline(0,1)
    }
}

#' Compute Bootstrap-Summary for 'MCResultBCa' Objects.
#' 
#' Function computes the bootstrap summary for objects of class 'MCResultBCa'.
#' 
#' @param .Object object of class 'MCResultBCa'
#' 
#' @return matrix of bootstrap results

MCResultBCa.bootstrapSummary<-function(.Object)
{
    bsum <- matrix(NA, nrow=2, ncol=4)
    colnames(bsum) <- c("global.est","bootstrap.mean","bias","bootstrap.se")
    rownames(bsum) <- c("Intercept","Slope")
    bsum[,"global.est"] <- .Object@glob.coef
    bsum[,"bootstrap.mean"] <- c(mean(.Object@B0), mean(.Object@B1))
    bsum["Intercept","bias"] <- (1/.Object@nsamples)*sum(.Object@B0-.Object@glob.coef[1])   #????
    bsum["Slope","bias"] <- (1/.Object@nsamples)*sum(.Object@B1-.Object@glob.coef[2])   #????
    bsum[,"bootstrap.se"] <- c(sd(.Object@B0), sd(.Object@B1))
    return(round(bsum, 5))
}

#' Caluculate Response 
#' 
#' Calculate predicted values for given values of the reference-method.
#' @param .Object object of class 'MCResultBCa'
#' @param x.levels numeric vector specifying values of the reference method for which prediction should be made
#' @param alpha significance level for confidence intervals
#' @param bootcimeth character string specifying the method to be used for bootstrap confidence intervals

MCResultBCa.calcResponse<-function(.Object, x.levels,alpha=0.05, bootcimeth=.Object@bootcimeth)
{
    stopifnot(is.numeric(alpha))
    stopifnot(alpha > 0 & alpha < 1)
    stopifnot(length(alpha) > 0)
    stopifnot(!is.na(x.levels))
    stopifnot(is.numeric(x.levels))
    stopifnot(length(x.levels) > 0)
    stopifnot(bootcimeth %in% c("Student","tBoot","quantile","BCa"))
    
    if (.Object@regmeth %in% c("PaBa", "PaBaLarge", "WDeming") & .Object@cimeth == "bootstrap" & bootcimeth=="tBoot")  
    stop(paste("It is impossible  to calculate the tBoot confidence bounds for ", .Object@regmeth, ".\n Please choose nested bootstrap.", sep=""))    

    npoints<-length(.Object@data[,"x"])
    nsamples <- .Object@nsamples
    cimeth <- .Object@cimeth

    b0 <- .Object@glob.coef[1]
    b1 <- .Object@glob.coef[2]

    sigmab0 <- .Object@glob.sigma[1]
    sigmab1 <- .Object@glob.sigma[2]

    B0jack <- .Object@B0jack
    B1jack <- .Object@B1jack

    B0 <- .Object@B0
    B1 <- .Object@B1

    sigmaB0 <- .Object@sigmaB0
    sigmaB1 <- .Object@sigmaB1

    xw <- .Object@xmean
    MX <- .Object@MX

    ylevels <- matrix(nrow=nsamples, ncol=length(x.levels))   # ylevels[k,j]= b0(k)+b1(k)*x.level[j]
    colnames(ylevels) <- paste("X", 1:length(x.levels), sep="")
    for(i in seq(along = x.levels)) ylevels[,i] <- B0+B1*x.levels[i]

    mresp <- matrix(ncol=5,nrow=length(x.levels))
    colnames(mresp) <- c("X","Y","Y.SE","Y.LCI","Y.UCI")
    mresp[,"X"] <- x.levels

    if(bootcimeth == "Student")
    {
        mresp[,"Y"] <- b0+b1*mresp[,"X"]
        mresp[,"Y.SE"] <- apply(ylevels, 2, sd, na.rm=TRUE)
        mresp[,"Y.LCI"] <- mresp[,"Y"]-qt(1-alpha/2, npoints-2)*mresp[,"Y.SE"]
        mresp[,"Y.UCI"] <- mresp[,"Y"]+qt(1-alpha/2, npoints-2)*mresp[,"Y.SE"]
    } 

    if(bootcimeth == "tBoot")
    {
        mresp[,"Y"]<-b0+b1*mresp[,"X"]
        if(cimeth == "nestedbootstrap")
        {
            mresp[,"Y.SE"] <- apply(ylevels, 2, sd, na.rm=TRUE)
        }
        else
        {
            mresp[,"Y.SE"] <- sqrt( .Object@para["Intercept", "SE"]^2 +
                                    .Object@para["Slope", "SE"]^2 * mresp[,"X"] * (mresp[,"X"]
                                    - 2 * xw))
        }
        Tstar <- mc.calcTstar(.Object, x.levels)
        tstar.low<-apply(Tstar,2,mc.calc.quant,alpha=alpha/2)
        tstar.up<-apply(Tstar,2,mc.calc.quant,alpha=1-alpha/2)
        mresp[,c("Y.LCI","Y.UCI")] <- c(mresp[,"Y"]-tstar.up*mresp[,"Y.SE"],mresp[,"Y"]-tstar.low*mresp[,"Y.SE"])
    } 

    if(bootcimeth == "BCa")
    {
        mresp[,"Y"] <- b0+b1*mresp[,"X"]
        mresp[,"Y.SE"]<-  NA
        ylevelsJack <- matrix(nrow=npoints, ncol=length(x.levels))   # ylevels[k,j]= b0(k)+b1(k)*x.level[j]
        colnames(ylevels) <- paste("X", 1:length(x.levels), sep="")
        for(j in seq_along(x.levels)) 
            ylevelsJack[,j] <- B0jack+B1jack*x.levels[j]

        for(k in seq(along=x.levels))
        {
            BCA <- mc.calc.bca(Xboot=ylevels[,k], Xjack=ylevelsJack[,k], xhat=mresp[k,"Y"], alpha=alpha)
            mresp[k,c("Y.LCI","Y.UCI")] <- BCA$CI
       }
    } 

    if(bootcimeth == "quantile")
    {
        mresp[,"Y"] <- b0+b1*mresp[,"X"]
        mresp[,"Y.LCI"] <- apply(ylevels,2,quantile, probs=alpha/2, type=3)
        mresp[,"Y.UCI"] <- apply(ylevels,2,quantile, probs=1-alpha/2, type=3)
    }
    
    return(mresp)
}

#' Print Regression-Analysis Summary for Objects of class 'MCResultBCa'.
#' 
#' Functions prints a summary of the regression-analysis for objects of class 'MCResultBCa'.
#' 
#' @param .Object object of class 'MCResultBCa'

MCResultBCa.printSummary<-function(.Object)
{
    regmeth <- .Object@regmeth
    cimeth <- .Object@cimeth
    bootcimeth <- .Object@bootcimeth
    regtext <- ""

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
    cat(paste("The confidence intervals are calculated with\n", cimeth," (", bootcimeth,") method.\n"), sep="")
    cat(paste("Confidence level: ", (1-.Object@alpha)*100, "%\n", sep=""))
    if(regmeth %in% c("Deming","WDeming")) 
        cat("Error ratio:", .Object@error.ratio)
    cat("\n\n------------------------------------------\n\n")
    cat(toupper(paste(regtext, "Fit:\n\n")))
    print(getCoefficients(.Object))
    cat("\n\n------------------------------------------\n\n")
    cat("BOOTSTRAP SUMMARY\n\n")
    print(bootstrapSummary(.Object))
    if(!is.na(.Object@rng.seed)) {
        cat("\nBootstrap results generated with fixed RNG settings.\n")
        cat(paste("RNG kind: ",.Object@rng.kind,"\nRNG seed: ",.Object@rng.seed,"\n",sep=""))
    }
    else
    {
        cat("\nBootstrap results generated with environment RNG settings.\n")
    }
    return(NULL)
}
