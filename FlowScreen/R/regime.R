#' Plot flow regime
#' 
#' This function plots the min, max, mean, and two user-defined quantiles of 
#' daily streamflow to provide visual summary of the flow regime. Area between 
#' the upper and lower quantile is shaded grey, the dark blue line represents
#' the mean daily discharge, gray line represents the median daily discharge,
#'  and the period of record daily maximum and minimum are
#' shown with the blue points.
#' @param TS output from \code{\link{create.ts}} containing a data.frame of flow
#'   time series
#' @param q Numeric vector of the upper and lower quantile values.  Default
#'   is c(0.9, 0.1).
#' @param by Character string indicating whether to plot the regime by day of 
#'   the hydrologic year (defined using \code{\link{create.ts}}) or by day of 
#'   the calendar year. Options are "doy" (calendar year) or "hdoy" (default, hydrologic year).
#' @param text optional character string for margin text, e.g. for station name, 
#'   location, or other notes. Set to NULL if no margin text is wanted, or set to "d" 
#'   to use default text containing the station ID, station name, and province/state 
#'   returned from \code{\link{station.info}}. 
#' @author Jennifer Dierauer
#' @export
#' @examples
#' # plot the flow regime of the Caniapiscau River
#' data(cania.sub.ts)
#' regime(cania.sub.ts)


regime <- function(TS, q=c(0.9, 0.1), text="d", by="hdoy") {
    
    opar <- graphics::par(no.readonly = TRUE)
    hyrstart <- as.numeric(subset(TS, TS$hmonth==1)$month[1])
    
    if (by=="hdoy") {
        doy <- as.factor(TS$hdoy)
        
    } else {
        doy <- as.factor(TS$doy)
    } 
    
    ### initialize array to be filled
    Qdoy <- array(data=NA, c(max(as.numeric(doy)),6))
    colnames(Qdoy)<- c("MaxQ", "MinQ", "MeanQ", "Q90", "Q10", "Median")
    
    ### calculate stats
    Qdoy[,1]<-tapply(TS$Flow, doy, max, na.rm=TRUE)
    Qdoy[,2]<-tapply(TS$Flow, doy, min, na.rm=TRUE)
    Qdoy[,3]<-tapply(TS$Flow, doy, mean, na.rm=TRUE)
    Qdoy[,4]<-tapply(TS$Flow, doy, stats::quantile, q[1], na.rm=TRUE)
    Qdoy[,5]<-tapply(TS$Flow, doy, stats::quantile, q[2], na.rm=TRUE)
    Qdoy[,6]<-tapply(TS$Flow, doy, stats::median, na.rm=T)
    
    ### set up polygon for inter-quantile shading
    mdoy <- as.numeric(unique(doy))
    xx <- c(1:max(as.numeric(doy)),max(as.numeric(doy)):1)
    yy <- c(Qdoy[,4],Qdoy[max(as.numeric(doy)):1,5])
    
    ### create plot
    graphics::par(mar=c(4,4,2,2))
    if (!is.null(text)) {graphics::par(oma=c(0,0,1,0))}
    
    yl1=expression(paste("Discharge (m" ^{3}, "/s)"))
    graphics::plot(Qdoy[,1], col="#6BAED6", type="p", pch=19, cex=0.5, xlab="", ylab="",
         xaxt="n")#max
    graphics::title(ylab=yl1, line=2)
    graphics::points(Qdoy[,2], col="#6BAED6", type="p", pch=19, cex=0.5) #min
    graphics::polygon(xx, yy, col="gray", border="#3182BD")
    graphics::points(Qdoy[,3],col="#08519C",type="l",lwd=2) #mean
    graphics::points(Qdoy[,6], col="gray50", type="l", lwd=2) #median
    
    axis_doy.internal(hyrstart)
    
    SYMnames <- c("maximum", paste(as.character(max(q)), "quantile"), "mean", "median", 
                  paste(as.character(min(q)), "quantile"), "minimum")
    SYMcol <- c("#6BAED6", "#3182BD", "#08519C", "gray50", "#3182BD", "#6BAED6")
    
    graphics::legend("topleft", legend = SYMnames, lwd = c(NA, 1, 2, 2, 1, NA), 
           lty = c(NA, 1, 1, 1, 1, NA),
           pch = c(19, NA, NA, NA, NA, 19), pt.cex = 0.5, col = SYMcol, cex=0.7)
    
    if (!is.null(text)) {
        
        if (text == "d") {
            stinfo <- station.info(TS, Plot=F)
            text <- paste("ID: ", stinfo[1],", NAME: ", stinfo[2],
                          ", PROV/STATE: ", stinfo[3], sep = "")
        }
        
        graphics::mtext(text, side=3, line=0, outer=T, cex=0.7)
    }
    
    on.exit(suppressWarnings(graphics::par(opar)))
}
