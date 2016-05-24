################################################################################
##
## mcrCompareFit.r
##
## Grafical comparison of estimates of two or more regression fits.
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
################################################################################

#' Graphical Comparison of Regression Parameters and Associated Confidence Intervals
#' 
#' Graphical comparison of regression parameters (intercept and slope) and their associated
#' 100(1-alpha)\% confidence intervals for multiple fitted models of 'MCResult' sub-classes.
#'
#' @param ... list of fitted models, i.e. objects of "MCResult" sub-classes.
#' @examples
#'      library("mcr")
#'      data("creatinine", package="mcr")
#'      fit.lr <- mcreg(as.matrix(creatinine), method.reg="LinReg", na.rm=TRUE)
#'      fit.wlr <- mcreg(as.matrix(creatinine), method.reg="WLinReg", na.rm=TRUE)
#'      compareFit( fit.lr, fit.wlr )
compareFit <- function( ... )
{
    models <- list(...)
    
    for(i in seq_along(models)) 
        stopifnot(class(models[[i]])[1] %in% c("MCResultAnalytical", "MCResultJackknife", "MCResultResampling", "MCResultBCa"))

    nm <- length(models)

    # For each model collect informations about slope, intercept

    slopes <- matrix(nrow=length(models),ncol=3)
    colnames(slopes)<-c("EST","LCI","UCI")

    intercepts<- matrix(nrow=length(models),ncol=3)
    colnames(intercepts)<-c("EST","LCI","UCI")

    # For each model collect informations about regression type and inference type
    axis.lab <- c()

    for(i in seq_along(models))
    {
        slopes[i,]<-getCoefficients(models[[i]])["Slope",c("EST","LCI","UCI")]
        intercepts[i,]<-getCoefficients(models[[i]])["Intercept",c("EST","LCI","UCI")]
        if(models[[i]]@regmeth %in% c("Deming","WDeming"))
        {
            if (class(models[[i]])[1]!="MCResultResampling") 
                axis.lab[i] <- paste(models[[i]]@regmeth," fit \n",models[[i]]@cimeth," CI\nlambda=",models[[i]]@error.ratio,sep="")
            else
                axis.lab[i] <- paste(   models[[i]]@regmeth," fit\n",models[[i]]@cimeth,"(",models[[i]]@bootcimeth,") CI\nlambda=",
                                        models[[i]]@error.ratio,sep="")
        } 
        else if(class(models[[i]])[1]!="MCResultResampling") 
            axis.lab[i] <- paste(models[[i]]@regmeth," fit\n",models[[i]]@cimeth," CI",sep="")
        else 
            axis.lab[i] <- paste(models[[i]]@regmeth," fit\n",models[[i]]@cimeth,"(",models[[i]]@bootcimeth,") CI",sep="")
    }# end of loop

    YLIM <- c(0.6,nm+0.3) # limits for grafic
    
    old.par <- par(c("mfrow", "oma", "mar"))

    par(mfrow=c(1,3),oma=c(2,2,2,2))

    par(mar=c(4, 0, 3, -0.1) + 0.1)
    
    plot(c(0,1),c(0,1),ylim=YLIM, xlab="",ylab="", col="white",yaxt="n",xaxt="n",bty="n",main="")
    mtext(side=2, line=0,at=(1:nm),axis.lab,adj=0,cex=0.8,las=1)

    par(mar=c(4, -0.1, 3, 2) + 0.1)

    delta <- as.numeric(colMeans(intercepts)[3]-colMeans(intercepts)[2])
    XLIM<-range(intercepts)+c(-1/2,1/2)*delta

    plot(0,0,xlim=XLIM,ylim=YLIM, xlab="",ylab="", col="white",yaxt="n",bty="n",main="Intercept",cex.axis=1)
    abline(h=1:nm, col="lightgrey",lty=2)
    points(intercepts[,"EST"],1:nm,pch=16)
    for(s in 1:nm) 
        arrows(intercepts[s, "LCI"],s,intercepts[s, "UCI"],s,angle=90,code=3,length=0.06)
    
    abline(v=0, col="red", lty=2)

    par(mar=c(4, 2, 3, 2) + 0.1)

    delta <- as.numeric(colMeans(slopes)[3]-colMeans(slopes)[2])
    XLIM <- range(slopes)+c(-1/2,1/2)*delta

    plot(0,0,xlim=XLIM,ylim=YLIM, xlab="",ylab="", col="white",yaxt="n",bty="n", main="Slope",cex.axis=1)
    abline(h=1:nm, col="lightgrey",lty=2)
    points(slopes[,"EST"],1:nm,pch=16)
    for (s in 1:nm) arrows(slopes[s, "LCI"],s,slopes[s, "UCI"],s,angle=90,code=3,length=0.06)    # code=2 - zweiseitig
        abline(v=1, col="red", lty=2)
        
    par(old.par)
} #end of function


