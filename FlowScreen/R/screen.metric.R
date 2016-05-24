#' Plot a metric with trend and change points
#' 
#' This function plots a time series of a streamflow metric with the prewhitened
#' linear trend and any detected changepoints in mean and variance. 
#' @param y Numeric vector with "times" attribute
#' @param ylabel Character string for the y-axis label
#' @param text optional character string for margin text, e.g. for station name, 
#'   location, or other notes.
#' @details This function plots detected changepoints as a vertical dashed line.
#'   The means on either side of a changepoint are plotted as solid black lines.
#'   If the temporal trend is significant (p-value < 0.1), the trend is plotted as 
#'   a blue or red line for an increasing or decreasing trend, respectively.
#'   The upper and lower 95% confidence bounds for the trend are represented by the
#'   dotted red or blue lines.  If a trend is not significant, it is not plotted.
#' @return Returns a list containing results from the trend and changepoint
#'   analysis. This list has the following elements:
#'   \itemize{
#'     \item slope - Numeric vector containing the intercept and slope of the 
#'       prewhitened linear trend computed with \code{\link[zyp]{zyp.trend.vector}}
#'       using Yue Pilon's method
#'     \item ci1 - numeric vector containing the intercept and slope of the
#'       upper confidence bound. See \code{\link[zyp]{confint.zyp}}
#'     \item ci2 - numeric vector of length 2 containing the intercept and slope
#'       of the lower confidence bound. See \code{\link[zyp]{confint.zyp}}
#'     \item pval - numeric value indicatng the significance value of the detected
#'       trend, Kendall test computed within \code{\link[zyp]{zyp.trend.vector}}
#'     \item cpts - numeric vector of changepoints if any are found, computed 
#'       with \code{\link[changepoint]{cpt.meanvar}}
#'     \item means - numeric vector of means computed with 
#'       \code{\link[changepoint]{cpt.meanvar}}
#'   }
#' @author Jennifer Dierauer
#' @seealso See \code{\link{screen.summary}} to create a summary screening plot of 
#'   high flow, low flow, or baseflow metrics.
#'   
#'   See \code{\link{metrics.all}} to calculate 30 different streamflow metrics at once.
#'   The \code{\link{screen.metric}} function could then be used to loop through the metrics and 
#'   create an individual plot for each.
#' @export
#' @examples
#' data(cania.sub.ts)
#' 
#' # calculate and plot the annual maximum series
#' res <- pk.max(cania.sub.ts)
#' res1 <- screen.metric(res, ylabel="Q (m3/s)", 
#' text="Caniapiscau River, Annual Maximum Series")
#' 
#' # calculate and plot the annual minimum series
#' res <- MAMn(cania.sub.ts, n=1)
#' res1 <- screen.metric(res, ylabel="Discharge (m3/s)", 
#' text="Caniapiscau River, Annual Minimum Series")

screen.metric <- function(y, ylabel="", text=NULL) {
    
    opar <- graphics::par(no.readonly = TRUE)
    
    MyY <- y
    MyX <- attr(MyY, "times")
    
    ### set y axis limits
    MyYlims <- c(0, ceiling(1.2*max(MyY, na.rm=TRUE)))
    
    ### format x values to work with plotting of sen slopes and change points
    if(length(MyX[1]) > 5) {
        Year <- as.numeric(format(MyX, "%Y"))
        Year1 <- min(Year)
        YearEnd <- max(Year)
        Start <- as.Date(paste(Year1, "-01-01", sep=""))
        MyX.mod <- c(1:length(MyX))
        for (j in 1:length(MyX)) {MyX.mod[j] <- (MyX[j]-Start)}
    } else {
        MyX.mod <- c(1:length(MyX))
        for (j in 2:length(MyX)) {MyX.mod[j] <- (as.numeric(MyX[j])-as.numeric(MyX[1]) + 1)}
    }

    #plot time series
    if (!is.null(text)) {graphics::par(oma=c(0,0,1,0))} else {graphics::par(oma=c(0,0,0,0))}
    
    graphics::par(mar=c(4,4,2,2))
    graphics::plot(MyX.mod, MyY, ylab=ylabel, xaxt="n", xlab="", type="p",
         pch=19, ylim=MyYlims)
    
    # add optional margin text
    if (!is.null(text)) {graphics::mtext(text, side=3, line=0, outer=T, cex=0.7)}
    
    #add x-axis ticks and labels
    if(nchar(as.character(MyX[1]), type="chars") > 4) {
        Years <- substr(MyX, 1, 4)
        Year1 <- as.numeric(Years[1])
        YearEnd <- as.numeric(Years[length(Years)])
        myticks <- seq(from=0, to=365*(YearEnd-Year1), by=365)
        graphics::axis(1, at=myticks, labels=c(Year1:YearEnd))
    } else {graphics::axis(1, at=1:((max(MyX.mod))), labels=c(min(MyX):max(MyX)))}
    
    # used zyp.sen and confint.zyp to get the intercept and slope of the 
    # upper and lower confidence intervals
    mrange <- max(MyY, na.rm=T) - min(MyY, na.rm=T)
    if (mrange > 0) {
        slope <- zyp::zyp.sen(MyY~MyX.mod)
        ci <- zyp::confint.zyp(slope)
        ci1 <- ci[,2]
        ci2 <- ci[,1]
        
        # trend intercept, slope, and p-value are calculated as the prewhitened
        # linear trend with the Yue Pilon method.
        res <- zyp::zyp.trend.vector(MyY, x=MyX.mod, method="yuepilon")
        slope <- c(res[[11]], res[[2]])
        pval <- res[[6]]
        
    } else {
        slope <- NA
        ci1 <- NA
        ci2 <- NA
        pval <- NA
    }

    #add trend line if p-value is less than 0.1
    if (!is.na(pval) && pval <= 0.1) {
        mcol <- ifelse(slope[2] < 0, "darkred", "darkblue")
        mlwd <- ifelse(pval <= 0.05, ifelse(pval<=0.01, 3, 2), 1)
        graphics::abline(coef=slope, col=mcol, lwd=mlwd)
        graphics::abline(coef=ci1, lty=3, col=mcol)
        graphics::abline(coef=ci2, lty=3, col=mcol)
        
        ypos <- ifelse(mean(MyY[1:3]) <= 0.5*MyYlims[2], 0.95*MyYlims[2], 1.1*MyYlims[1])
 
        mypvalue <- round(as.numeric(pval), digits=2)
        
        if (mypvalue == 0) {
            graphics::text(MyX.mod[1], ypos, paste("Trend p-value < 0.01"),col=mcol, 
                 adj=c(0,0))
        } else {
            graphics::text(MyX.mod[1], ypos, paste("Trend p-value =", mypvalue),col=mcol, 
                adj=c(0,0))
        }
    }
    
    # find change points
    if (length(MyY) > 3) {
        out <- suppressWarnings(changepoint::cpt.meanvar(as.numeric(MyY), 
                                            penalty="Asymptotic",
                                            pen.value=0.05,
                                            method="BinSeg"))
        MyCpts <- out@cpts
        MyMeans <- out@param.est$mean
        
        #add vertical lines for change points and horizontal lines for means
        NumPoints <- length(MyX.mod)
        MyCpts <- MyCpts[MyCpts > 3 & MyCpts < (NumPoints-3)] ## remove cpts at end and beginning
        
        if (length(MyCpts) > 0){
            
            for (j in 1:length(MyCpts)) {
                graphics::abline(v=MyX.mod[MyCpts[j]], lwd=2, lty=5)
                ifelse(j==1, xpts <- c(-10, MyX.mod[MyCpts[j]]),
                       xpts <- c(MyX.mod[MyCpts[j-1]], MyX.mod[MyCpts[j]]))
                ypts <- c(MyMeans[j], MyMeans[j])
                graphics::points(xpts, ypts, type="l", lwd=2)
            }
            
            xpts <- c(xpts[2], 1.1*max(MyX.mod))
            ypts <- c(MyMeans[length(MyMeans)], MyMeans[length(MyMeans)])
            graphics::points(xpts, ypts, type="l", lwd=2)
        }
    }
    
    #compile results for return
    out <- list(Slope=slope, ci1=ci1, ci2=ci2, pval=pval, cpts=MyCpts, means=MyMeans)
    return(out)
    
    on.exit(suppressWarnings(graphics::par(opar)))
}

