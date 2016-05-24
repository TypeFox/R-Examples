cdf.plot <- function(cdfest, units.cdf="Percent", type.cdf="Continuous",
   logx="", xlbl=NULL, ylbl="Percent", ylbl.r=NULL, figlab=NULL, legloc="BR",
   confcut=5, conflev=95, cex.main=1.2, ...) {

################################################################################
# Function: cdf.plot
# Purpose: Plot the CDF and associated confidence limits
# Programmers: Tony Olsen
#              Tom Kincaid
# Date: March 26, 2007
# Last Revised: June 23, 2014
# Description:
#   This function creates a CDF plot.  Input data for the plots is provided by a
#   data frame utilizing the same structure as the data frame named "CDF" that
#   is included in the output object produced by function cont.analysis, but the
#   data frame includes only the values for a single CDF.  Confidence limits for
#   the CDF also are plotted.
# Arguments:
#   cdfest = data frame utilizing the same structure as the data frame named
#     "CDF" that is included in the output object produced by function
#     cont.analysis.  The data frame must contain only a single cdf estimate.
#   units.cdf = indicator for the type of units in which the CDF is plotted,
#     where "Percent" means the plot is in terms of percent of the population,
#     and "Units" means the plot is in terms of units of the population.  The
#     default is "Percent".
#   type.cdf = character string consisting of the value "Continuous" or
#     "Ordinal" that controls the type of CDF plot for each indicator.  The
#     default is "Continuous".
#   logx = character string consisting of the value "" or "x" that controls
#     whether the x axis uses the original scale ("") or the base 10 logarithmic
#     scale ("x").  The default is "".
#   xlbl = character string providing the x-axis label.  If this argument equals
#     NULL, then the indicator name is used as the label.  The default is NULL.
#   ylbl =  character string providing the y-axis label.  The default is
#     "Percent".
#   ylbl.r = character string providing the label for the right side y-axis,
#     where NULL means a label is not created, and "Same" means the label is the
#     same as the left side label (i.e., argument ylbl).  The default is NULL.
#   figlab = character string providing the plot title.  The default is NULL.
#   legloc = indicator for location of the plot legend, where "BR" means bottom
#    right, "BL" means bottom left, "TR" means top right, and "TL" means top
#    left.  The default is "BR". 
#   confcut = numeric value that controls plotting confidence limits at the CDF
#     extremes.  Confidence limits for CDF values (percent scale) less than
#     confcut or greater than 100 minus confcut are not plotted.  A value of
#     zero means confidence limits are plotted for the complete range of the
#     CDF.  The default is 5.
#   conflev = numeric value of the confidence level used for confidence limits.
#     The default is 95.
#   cex.main = expansion factor for the plot title.  The default is 1.2.
#   ... = additional arguments passed to the plot function.
# Output:
#   A plot of the CDF and its associated confidence limits.
# Other Functions Required:
#   interp.cdf - interpolate CDF values at a set of percentiles
#   interp.axis - create right side y-axis labels for a CDF plot
# Example:
#   mysiteID <- paste("Site", 1:100, sep="")
#   mysites <- data.frame(siteID=mysiteID, Active=rep(TRUE, 100))
#   mysubpop <- data.frame(siteID=mysiteID, All.Sites=rep("All Sites",100),
#      Resource.Class=rep(c("Good","Poor"), c(55,45)))
#   mydesign <- data.frame(siteID=mysiteID, wgt=runif(100, 10, 100),
#      xcoord=runif(100), ycoord=runif(100), stratum=rep(c("Stratum1",
#      "Stratum2"), 50))
#   ContVar <- rnorm(100, 10, 1)
#   mydata.cont <- data.frame(siteID=mysiteID, ContVar=ContVar)
#   mypopsize <- list(All.Sites=c(Stratum1=3500, Stratum2=2000),
#      Resource.Class=list(Good=c(Stratum1=2500, Stratum2=1500),
#      Poor=c(Stratum1=1000, Stratum2=500)))
#   myanalysis <- cont.analysis(sites=mysites, subpop=mysubpop, design=mydesign,
#      data.cont=mydata.cont, popsize=mypopsize)
#   keep <- myanalysis$CDF$Type == "Resource.Class" & 
#      myanalysis$CDF$Subpopulation == "Good"       
#   par(mfrow=c(2,1))       
#   cdf.plot(myanalysis$CDF[keep,], xlbl="ContVar",
#      ylbl="Percent of Stream Length", ylbl.r="Stream Length (km)",
#      figlab="Estimates for Resource Class: Good")
#   cdf.plot(myanalysis$CDF[keep,], xlbl="ContVar",
#      ylbl="Percent of Stream Length", ylbl.r="Same)",
#      figlab="Estimates for Resource Class: Good")
################################################################################

# Set graphical parameter values

op <- par(mgp=c(1.7,0.6,0), mar=c(3,3,2,4)+0.1)

# Create the data frame of values to be plotted and set the y-axis limits
# The data frame structure follows: column 1: x-axis values
#                                   column 2: CDF estimates for left y-axis
#                                   column 3: lower confidence limit values
#                                   column 4: upper confidence limit values
#                                   column 5: CDF estimates for right y-axis

if(units.cdf == "Percent") {
   cdfdata <- cdfest[,c(4,6,8,9,10)]
} else if(units.cdf == "Units"){
   cdfdata <- cdfest[,c(4,10,12,13,6)]
} else{
   stop(paste("\nThe choice of units for the CDF must be either \"Percent\" or \"Units\". The value \nsupplied for argument units.cdf was: \"", units.cdf, "\".\n", sep=""))
}

# Restrict confidence limits to lie between confcut and 100-confcut percent

pctval <- c(confcut, 100-confcut)
tvalue <- cdfest[,6] >= pctval[1] & cdfest[,6] <= pctval[2]
x <-  interp.cdf(pctval, cdfest[,6], cdfdata[,1])
ylow <- interp.cdf(pctval, cdfest[,6], cdfdata[,3])
yhi <- interp.cdf(pctval, cdfest[,6], cdfdata[,4])

# Set the left side y-axis limits

if(units.cdf == "Percent") {
  ylimit <- c(0,100)
} else if(units.cdf == "Units"){
  ylimit <- pretty(c(min(c(cdfdata[,2], ylow)), max(c(cdfdata[,2], yhi))))
  ylimit <- ylimit[c(1, length(ylimit))]
}

# Plot the CDF for a continuous indicator

if(type.cdf == "Continuous") {
   plot(cdfdata[,1], cdfdata[,2], type="l", ylim=ylimit, xlab=xlbl, ylab=ylbl,
        log=logx, ...)

   # Plot confidence limits

   value <- c(x[1], cdfdata[,1][tvalue], x[2])
   lower <- c(ylow[1], cdfdata[,3][tvalue], ylow[2])
   upper <- c(yhi[1], cdfdata[,4][tvalue], yhi[2])
   lines(value, lower, lty=3, lwd=1.5)
   lines(value, upper, lty=3, lwd=1.5)

# Plot the CDF for an ordinal (count) indicator

} else if(type.cdf == "Ordinal") {
   x <- rep(cdfdata[,1], each=2)[-1]
   y <- rep(cdfdata[,2], each=2)
   tmp <- cbind(matrix(c(x,x[length(x)]),ncol=2,byrow=TRUE),rep(NA,nrow(cdfdata)))
   x <- as.vector(t(tmp))
   tmp <- cbind(matrix(y,ncol=2,byrow=TRUE),rep(NA,nrow(cdfdata)))
   y <- as.vector(t(tmp))
   plot(x, y, type="l", ylim=ylimit, xlab=xlbl, ylab=ylbl, ...)

   # Plot confidence limits

   len <- length(cdfdata[,1][tvalue])
   if(len > 1) {
      value <- rep(cdfdata[,1][tvalue], each=2)[-1]
      tmp <- cbind(matrix(c(value,value[length(value)]),ncol=2,byrow=TRUE),
                   rep(NA,len))
      value <- as.vector(t(tmp))
      len <- length(cdfdata[,4][tvalue])
      if(len > 1) {
         lower <- rep(cdfdata[,3][tvalue], each=2)
         tmp <- cbind(matrix(lower,ncol=2,byrow=TRUE),rep(NA,len))
         lower <- as.vector(t(tmp))
         upper <- rep(cdfdata[,4][tvalue], each=2)
         tmp <- cbind(matrix(upper,ncol=2,byrow=TRUE),rep(NA,len))
         upper <- as.vector(t(tmp))
         lines(value,lower,lty=3, lwd=1.5)
         lines(value,upper,lty=3, lwd=1.5)
      }
   }
} else {
   stop(paste("\nThe type of CDF must be either \"Continuous\" or \"Ordinal\". The value supplied \nfor argument type.cdf was: \"", type.cdf, "\".\n", sep=""))
}

# Create the plot title 

title(figlab, line=1, cex.main=cex.main)

# Create the plot legend

rx <- range(par("usr")[1:2], cdfdata[,1])
ry <- range(par("usr")[3:4], cdfdata[,2])
if(legloc=="BR") {
   xjust <- 1; yjust <- 0; 
   legx <- rx[2]; legy <- ry[1]
} else if(legloc=="BL") {
   xjust <- 0; yjust <- 0; 
   legx <- rx[1]; legy <- ry[1]
} else if(legloc=="TR") {
   xjust <- 1; yjust <- 1; 
   legx <- rx[2]; legy <- ry[2]
} else if(legloc=="TL") {
   xjust <- 0; yjust <- 1; 
   legx <- rx[1]; legy <- ry[2]
}
legend(x=legx, y=legy, xjust=xjust, yjust=yjust,
       legend=c("CDF Estimate",paste(conflev,"% Confidence Limits", sep="")),
       lty=c(1,3), lwd=c(1,1.5), bty="n", cex=1)

# If requested, create the right side y-axis labels

if(!is.null(ylbl.r)) {
   yl.lab <- seq(par("yaxp")[1], par("yaxp")[2], len=par("yaxp")[3]+1)
   if(ylbl.r == "Same") {
      axis(side=4, at=yl.lab, labels=yl.lab)
      mtext(ylbl, side=4, line=2, cex=par("cex"))
   } else {
      yr.lab <- interp.axis(yl.lab, cdfdata[,2], cdfdata[,5])
      axis(side=4, at=yl.lab, labels=as.character(round(yr.lab)))
      mtext(ylbl.r, side=4, line=2, cex=par("cex"))
   }
}

# Reset graphical parameter values

par(op)

invisible(NULL)
}
