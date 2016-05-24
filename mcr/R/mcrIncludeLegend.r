###############################################################################
##
## mcrIncludeLegend.r
##
## Draw legend for regression plot. 
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


#' Include Legend 
#' 
#' Include legend in regression plot (function \code{plot()}) or in bias plot (function \code{plotBias  ()}) with two or more lines.
#'
#' @param models list of length n with Objects of class "MCResult".
#' @param digits number of digits in Coefficients.
#' @param design type of legend design. There are two possible designs: "1" and "2" (See example).
#' @param place place for Legend: "topleft","topright","bottomleft" or "bottomright".
#' @param colors vector of length n with color of regression lines.
#' @param lty vector of length n with type of regression lines.
#' @param lwd  vector of length n with thickness of regression lines.
#' @param box.lty box line-type
#' @param cex numeric value representing the plotting symbol magnification factor
#' @param bg the background-color of the legend box
#' @param inset inset distance(s) from the margins as a fraction of the plot region when legend is placed by keyword.
#' @param bias logical value. If bias = TRUE, it will be drawn a legend for \code{plotBias()} function.
#' @param model.names legend names for different models. If NULL the regression type will be used.
#' @param ... other parameters of function legend().
#' @return Legend in plot.
#' @aliases includeLegend
#' @seealso \code{\link{plot.mcr}}, \code{\link{plotBias}}, \code{\link{plotResiduals}}, \code{\link{plotDifference}}, \code{\link{compareFit}}
#' @examples
#' #library("mcr")
#'
#'  data(creatinine,package="mcr")
#'  x <- creatinine$serum.crea
#'  y <- creatinine$plasma.crea
#'
#'  m1 <- mcreg(x,y,method.reg="Deming", mref.name="serum.crea",
#'                                         mtest.name="plasma.crea", na.rm=TRUE)
#'  m2 <- mcreg(x,y,method.reg="WDeming", method.ci="jackknife", 
#'                                          mref.name="serum.crea",
#'                                          mtest.name="plasma.crea", na.rm=TRUE)
#'
#'  plot(m1,  XLIM=c(0.5,3),YLIM=c(0.5,3), Legend=FALSE, 
#'                           Title="Deming vs. weighted Deming regression", 
#'                           Points.pch=19,ci.area=TRUE, ci.area.col=grey(0.9),
#'                           identity=FALSE, Grid=FALSE, Sub="")
#'  plot(m2, ci.area=FALSE, ci.border=TRUE, ci.border.col="red3", 
#'                           reg.col="red3", Legend=FALSE,add=TRUE, 
#'                           Points=FALSE, identity=FALSE, Grid=FALSE)
#'
#'  includeLegend(place="topleft",models=list(m1,m2), 
#'                           colors=c("darkblue","red"), design="1", digits=2)
includeLegend <- function(models = list(), digits = 2,design = paste(1:2),
                  place = c("topleft","topright","bottomleft","bottomright"),
                  colors, lty=rep(1, length(models)), lwd=rep(2, length(models)),
                  box.lty = "blank", cex=0.8, bg = "white",
                  inset = c(0.01,0.01),
                  bias=FALSE,model.names=NULL,...) {

        nm <- length(models)

        stopifnot(length(colors)==nm)
        stopifnot(length(lty)==nm)
        stopifnot(length(lwd)==nm)
        design <- match.arg(design)
        place <- match.arg(place)
        stopifnot(is.element(design, c("1","2")))
        stopifnot(is.element(place,c("topleft","topright","bottomleft","bottomright")))

        RegType <- c()
        RegEq <- c()
        CImethod <- c()
    
        
        for (i in seq_along(models)) {
            x <- models[[i]]
                
            if(is.null(model.names)) {
         	    if (x@regmeth == "LinReg") {
        		        RegType[i] <-"Linear Regression"
        	    } else if (x@regmeth == "WLinReg") {
        		    RegType[i] <-"Weighted Linear Regression"
        	    } else if (x@regmeth == "Deming") {
        		    RegType[i] <-"Deming Regression"
        	    } else if (x@regmeth == "WDeming") {
                    RegType[i] <-"Weighted Deming Regression"
                } else { 
                    RegType[i]  <- "Passing Bablok Regression"
                }
            } else RegType[i] <- model.names[i]
          
        RegEq[i] <- paste(round(x@para["Intercept","EST"],digits)," + ",
                          round(x@para["Slope","EST"],digits)," * ",x@mnames[1],sep="")
                      
        if (as.character(class(x)) %in% c("MCResultResampling","MCResultBCa" )) {
            CImethod[i] <- paste(x@cimeth," (",x@bootcimeth,")", sep="")
        } else {
            CImethod[i] <- x@cimeth   
        }
    }
        
    if (design == "1") {
        newlty <-rep(0, 2*nm)
        newlty[seq(1, 2*nm-1, by=2)] <-lty
        newlwd <-rep(0, 2*nm)
        newlwd[seq(1, 2*nm-1, by=2)] <-lwd
        newcolors <-rep("white", 2*nm)
        newcolors[seq(1, 2*nm-1, by=2)] <-colors
        newlegend <- c()
        newlegend[seq(1, 2*nm-1, by=2)] <- RegType
        if (bias == FALSE) {
            newlegend[seq(2, 2*nm, by=2)] <- RegEq
        } else {
            newlegend[seq(2, 2*nm, by=2)] <- CImethod
        }
    } else {
        newlty <- lty
        newlwd <- lwd
        newcolors <- colors
        if (bias == FALSE) {
            newlegend <- paste(RegType,": ",RegEq, sep="")
        } else {
            newlegend <- paste(RegType,": ",CImethod, sep="")
        }
    }            
                            
    legend(place, lty=newlty, lwd=newlwd, col=newcolors, legend=newlegend,cex=cex, bg=bg,inset=inset, box.lty=box.lty,...)
}
