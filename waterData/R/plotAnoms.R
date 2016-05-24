#' Function to plot hydrologic time series and anomalies
#'
#' @name plotAnoms
#' @title Plots streamflow anomalies
#' @param data is the anomaly list from the function 
#' \link{compAnom}.
#' @param ... further arguments to be passed to plotting method (see \link{par}).
#' @return a plot.
#' @keywords hplot ts multivariate
#' @export
#' @examples
#' q05054000.85 <- importDVs("05054000", sdate="1985-01-01", edate="2010-09-30")
#' anoms05054000 <- compAnom(q05054000.85, which=1)
#' plotAnoms(anoms05054000)
plotAnoms<- function(data, ...) {
  if ( ! is.data.frame(data) ) {
    if ( length(data) == 4 ) {
      par(mfrow=c(2,2), cex.lab=.9, las=1, tcl=0.5, xaxs="r", yaxs="r", 
          mar=c(4, 5, 0, 1) + 0.1,mgp=c(4, 1, 0),oma=c(0,0,1,0)+.1)
      plot(data[[1]]$dates, data[[1]]$val, log="y", type="l", xlab="",
           ylab="Streamflow, cubic feet per second", cex.lab=0.6, ...)
      plot(data[[1]]$dates, data[[1]]$ltfa, type="l", xlab="", 
           ylab="Long-term flow anomaly", cex.lab=0.6, ...)
      ltfa.sd <- round(sd(data[[1]]$ltfa, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", ltfa.sd, sep=""),side=3, line=-1, 
            cex=0.5)
      plot(data[[1]]$dates, data[[1]]$mtfa, type="l", xlab="", 
           ylab="Medium-term flow anomaly", cex.lab=0.6, ...)
      mtfa.sd<-round(sd(data[[1]]$mtfa, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", mtfa.sd, sep=""), side=3, line=-1, 
            cex=0.5)
      plot(data[[1]]$dates, data[[1]]$stfa, type="l", xlab="", 
           ylab="Short-term flow anomaly", cex.lab=0.6, ...)
      stfa.sd<-round(sd(data[[1]]$stfa, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", stfa.sd, sep=""), side=3, line=-1, 
            cex=0.5)
      plttxt<-paste("Site", as.character(data[[1]]$staid[1]), "LTFA is", data[[2]],
                    "days", "MTFA is", data[[3]], "days", "STFA is", data[[4]],
                    sep=" ")
      mtext(plttxt, outer=TRUE, side=3, cex=.8)
    }
    else if ( length(data) == 3 ) {
      par(mfrow=c(2,2), cex.lab=.9, las=1, tcl=0.5, xaxs="r", yaxs="r", 
      mar=c(4, 5, 0, 1) + 0.1,mgp=c(4, 1, 0),oma=c(0,0,1,0)+.1)
      plot(data[[1]]$dates, data[[1]]$val, log="y", type="l", xlab="",
           ylab="Streamflow, cubic feet per second", cex.lab=0.6, ...)
      plot(data[[1]]$dates, data[[1]]$mtfa, type="l", xlab="", 
           ylab="Medium-term flow anomaly", cex.lab=0.6, ...)
      mtfa.sd<-round(sd(data[[1]]$mtfa, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", mtfa.sd, sep=""), side=3, line=-1, 
            cex=0.5)
      plot(data[[1]]$dates, data[[1]]$stfa, type="l", xlab="", 
           ylab="Short-term flow anomaly", cex.lab=0.6, ...)
      stfa.sd<-round(sd(data[[1]]$stfa, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", stfa.sd, sep=""), side=3, line=-1, 
            cex=0.5)
      plttxt<-paste("Site", as.character(data[[1]]$staid[1]),
                    "MTFA is", data[[2]], "days", "STFA is", data[[3]],
                    sep=" ")
      mtext(plttxt, outer=TRUE, side=3, cex=.8)
    }
    else if ( length(data) == 6 ) {
      par(mfrow=c(3,2), cex.lab=.9, las=1, tcl=0.5, xaxs="r", yaxs="r", 
          mar=c(4, 5, 0, 1) + 0.1,mgp=c(4, 1, 0),oma=c(0,0,1,0)+.1)
      plot(data[[1]]$dates, data[[1]]$val, log="y", type="l", xlab="",
           ylab="Streamflow, cubic feet per second", cex.lab=0.6, ...)
      plot(data[[1]]$dates, data[[1]]$tenyranom, type="l", xlab="", 
           ylab="Ten-year flow anomaly", cex.lab=0.6, ...)
      tfa.sd <- round(sd(data[[1]]$tenyranom, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", tfa.sd, sep=""),side=3, line=-1, 
            cex=0.5)
      plot(data[[1]]$dates, data[[1]]$fiveyranom, type="l", xlab="", 
           ylab="Five-year flow anomaly", cex.lab=0.6, ...)
      ffa.sd <- round(sd(data[[1]]$fiveyranom, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", ffa.sd, sep=""),side=3, line=-1, 
            cex=0.5)
      plot(data[[1]]$dates, data[[1]]$annualanom, type="l", xlab="", 
         ylab="One-year flow anomaly", cex.lab=0.6, ...)
      ltfa.sd <- round(sd(data[[1]]$annualanom, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", ltfa.sd, sep=""),side=3, line=-1, 
            cex=0.5)
      plot(data[[1]]$dates, data[[1]]$seasanom, type="l", xlab="", 
           ylab="Seasonal flow anomaly", cex.lab=0.6, ...)
      mtfa.sd<-round(sd(data[[1]]$seasanom, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", mtfa.sd, sep=""), side=3, line=-1, 
            cex=0.5)
      plot(data[[1]]$dates, data[[1]]$dailyanom, type="l", xlab="", 
           ylab="Daily flow anomaly", cex.lab=0.6, ...)
      stfa.sd<-round(sd(data[[1]]$dailyanom, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", stfa.sd, sep=""), side=3, line=-1, 
            cex=0.5)
      plttxt<-paste("Site", as.character(data[[1]]$staid[1]), 
                    "Ten-year, five-year, one-year, seasonal, and daily anomalies",
                    sep=" ")
      mtext(plttxt, outer=TRUE, side=3, cex=.8)
    }
    else {
      stop("Data supplied must be a list of length 3, 4, or 6 generated by the 
           compQAnom function.")
    }
  }
  else if (is.data.frame(data)) {
    if ( length(data) == 6 ) {
      par(mfrow=c(2,2), cex.lab=.9, las=1, tcl=0.5, xaxs="r", yaxs="r", 
          mar=c(4, 5, 0, 1) + 0.1,mgp=c(4, 1, 0),oma=c(0,0,1,0)+.1)
      plot(data$dates, data$val, log="y", type="l", xlab="",
           ylab="Streamflow, cubic feet per second", cex.lab=0.6, ...)
      plot(data$dates, data$ltfa, type="l", xlab="", 
           ylab="Long-term flow anomaly", cex.lab=0.6, ...)
      ltfa.sd <- round(sd(data$ltfa, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", ltfa.sd, sep=""),side=3, line=-1, 
            cex=0.5)
      plot(data$dates, data$mtfa, type="l", xlab="", 
           ylab="Medium-term flow anomaly", cex.lab=0.6, ...)
      mtfa.sd<-round(sd(data$mtfa, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", mtfa.sd, sep=""), side=3, line=-1,
            cex=0.5)
      plot(data$dates, data$stfa, type="l", xlab="", 
           ylab="Short-term flow anomaly", cex.lab=0.6, ...)
      stfa.sd<-round(sd(data$stfa, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", stfa.sd, sep=""), side=3, line=-1,
            cex=0.5)
      plttxt<-paste("Site", as.character(data$staid[1]), sep=" ")
      mtext(plttxt, outer=TRUE, side=3, cex=.8)
    }
    else if ( length(data) == 5 ) {
      par(mfrow=c(2,2), cex.lab=.9, las=1, tcl=0.5, xaxs="r", yaxs="r", 
          mar=c(4, 5, 0, 1) + 0.1,mgp=c(4, 1, 0),oma=c(0,0,1,0)+.1)
      plot(data$dates, data$val, log="y", type="l", xlab="",
           ylab="Streamflow, cubic feet per second", cex.lab=0.6, ...)
      plot(data$dates, data$mtfa, type="l", xlab="", 
           ylab="Medium-term flow anomaly", cex.lab=0.6, ...)
      mtfa.sd<-round(sd(data$mtfa, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", mtfa.sd, sep=""), side=3, line=-1, 
            cex=0.5)
      plot(data$dates, data$stfa, type="l", xlab="", 
           ylab="Short-term flow anomaly", cex.lab=0.6, ...)
      stfa.sd<-round(sd(data$stfa, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", stfa.sd, sep=""), side=3, line=-1, 
            cex=0.5)
      plttxt<-paste("Site", as.character(data$staid[1]), sep=" ")
      mtext(plttxt, outer=TRUE, side=3, cex=.8)
    }
    else if ( length(data) == 8 ) {
      par(mfrow=c(3,2), cex.lab=.9, las=1, tcl=0.5, xaxs="r", yaxs="r", 
          mar=c(4, 5, 0, 1) + 0.1,mgp=c(4, 1, 0),oma=c(0,0,1,0)+.1)
      plot(data$dates, data$val, log="y", type="l", xlab="",
           ylab="Streamflow, cubic feet per second", cex.lab=0.6, ...)
      plot(data$dates, data$tenyranom, type="l", xlab="", 
           ylab="Ten-year flow anomaly", cex.lab=0.6, ...)
      tfa.sd <- round(sd(data$tenyranom, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", tfa.sd, sep=""),side=3, line=-1, 
            cex=0.5)
      plot(data$dates, data$fiveyranom, type="l", xlab="", 
           ylab="Five-year flow anomaly", cex.lab=0.6, ...)
      ffa.sd <- round(sd(data$fiveyranom, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", ffa.sd, sep=""),side=3, line=-1, 
            cex=0.5)
      plot(data$dates, data$annualanom, type="l", xlab="", 
           ylab="One-year flow anomaly", cex.lab=0.6, ...)
      ltfa.sd <- round(sd(data$annualanom, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", ltfa.sd, sep=""),side=3, line=-1, 
            cex=0.5)
      plot(data$dates, data$seasanom, type="l", xlab="", 
           ylab="Seasonal flow anomaly", cex.lab=0.6, ...)
      mtfa.sd<-round(sd(data$seasanom, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", mtfa.sd, sep=""), side=3, line=-1, 
            cex=0.5)
      plot(data$dates, data$dailyanom, type="l", xlab="", 
           ylab="Daily flow anomaly", cex.lab=0.6, ...)
      stfa.sd<-round(sd(data$dailyanom, na.rm=TRUE), digits=3)
      mtext(paste("Standard deviation is ", stfa.sd, sep=""), side=3, line=-1, 
            cex=0.5)
      plttxt<-paste("Site", as.character(data$staid[1]), 
                    "Ten-year, five-year, one-year, seasonal, and daily anomalies",
                  sep=" ")
      mtext(plttxt, outer=TRUE, side=3, cex=.8)
    }
    else {
      stop("Data supplied must be a data frame with 5, 6, or 8 columns.")
    }  
  }
  else {
    stop("Data must be a list or dataframe.")
  }
}