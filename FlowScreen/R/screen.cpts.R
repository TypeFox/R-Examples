#' Change point time series plot
#'
#' Compiles change point information for all metrics and outputs a daily flow
#' time series plot overlain with a bar plot of changepoint counts by year.
#' @param metrics output from \code{\link{metrics.all}}
#' @param type character indicating which type of metric to compile change points for.
#'   Options are "h" for high flow metrics, "l" for low flow metrics, "b" for baseflow 
#'   metrics, or "a" for all 30 metrics (10 high, 10 low, 10 baseflow).
#' @param text optional character string for margin text, e.g. for station name, 
#'   location, or other notes. Set to NULL if no margin text is wanted, or set to "d" 
#'   to use default text containing the station ID, station name, and province/state 
#'   returned from \code{\link{station.info}}.
#' @return When type="a", returns a data.frame of changepoint counts by metric 
#'   type and year.
#' @author Jennifer Dierauer
#' @seealso \code{\link{metrics.all}}
#' @export
#' @examples
#' # load results from metrics.all function for the Caniapiscau River
#' data(caniapiscau.res)
#' 
#' # plot changepoints for all metrics
#' screen.cpts(caniapiscau.res, type="l")

screen.cpts <- function(metrics, type="a", text=NULL) {
    
    opar <- graphics::par(no.readonly = TRUE)
    
    res <- metrics[[2]]
    TS <- metrics[[3]]
    cptsh <- list() 
    cptsl <- list()
    cptsb <- list()
    
    cmcolors <- grDevices::cm.colors(20, alpha=0.5)
    colh <- cmcolors[1]
    coll <- cmcolors[20]
    colb <- cmcolors[10]
    
    if (type == "h") {
        mtitle <- "Changepoints in High Flow Metrics"
        rylab <- "Number of Changepoints (max = 10)"
        types <- "Change Point Count"
        mcol <- colh
    }
    
    if (type == "l") {
        mtitle <- "Changepoints in Low Flow Metrics"
        rylab <- "Number of Changepoints (max = 10)"
        types <- "Change Point Count"
        mcol <- coll
    }
    
    if (type == "b") {
        mtitle <- "Changepoints in Baseflow Metrics"
        rylab <- "Number of Changepoints (max = 10)"
        types <- "Change Point Count"
        mcol <- colb
    }
    
    if (type == "a") {
        mtitle <- "Changepoints in All Metrics"
        rylab <- "Number of Changepoints (max = 30)"
        types <- c("Baseflow Changepoints", "Low Flow Changepoints", 
                   "High Flow Changepoints")
        mcol <- c(colb, coll, colh)
    }
    
    # put all change points in a list

    for (i in 1:10) {
        
        metname <- names(res[i])
        cptsh[[metname]] <- res[[i]]$cpts
        
    }
    for (i in 11:20) {
        
        metname <- names(res[i])
        cptsl[[metname]] <- res[[i]]$cpts
        
    }
    for (i in 21:30) {
        
        metname <- names(res[i])
        cptsb[[metname]] <- res[[i]]$cpts
        
    }
    
    ## identify start and end dates for plot, set axis limits
    Year1 <- TS$year[1]
    Start <- as.Date(paste(Year1, "-01-01", sep=""))
    YearLast <- TS$year[length(TS$year)]
    End <- as.Date(paste(YearLast, "-12-31", sep=""))
    myxlims <- c(Start, End)
    myylims <- c(0, 1.1 * max(TS$Flow))
    NumYears <- as.numeric(YearLast) - as.numeric(Year1) + 1
    
    Year_List <- seq(from=Year1, to=YearLast, by=1)
    
    CountsH <- c(numeric(0))
    CountsL <- c(numeric(0))
    CountsB <- c(numeric(0))
    
    for (i in 1:length(cptsh)) {
        cpts.sub <- cptsh[[i]]

        if (!is.na(cpts.sub[1])) {
            cpts.sub <- attr(cpts.sub, "times")
            for (j in 1:length(cpts.sub)) {
                ind <- as.numeric(substr(as.character(cpts.sub[j]), 1, 4))
                CountsH <- c(CountsH, ind)
            }      
        }
    }
    
    for (i in 1:length(cptsl)) {
        cpts.sub <- cptsl[[i]]
 
        if (!is.na(cpts.sub[1])) {
            cpts.sub <- attr(cpts.sub, "times")
            for (j in 1:length(cpts.sub)) {
                ind <- as.numeric(substr(as.character(cpts.sub[j]), 1, 4))
                CountsL <- c(CountsL, ind)
            }      
        }
    }
    
    for (i in 1:length(cptsb)) {
        cpts.sub <- cptsb[[i]]
 
        if (!is.na(cpts.sub[1])) {
            cpts.sub <- attr(cpts.sub, "times")
            
            for (j in 1:length(cpts.sub)) {
                ind <- as.numeric(substr(as.character(cpts.sub[j]), 1, 4))
                CountsB <- c(CountsB, ind)
            }      
        }
    }

    y1 <- expression(paste("Discharge (m" ^{3}, "/s)"))
    
    graphics::par(mar = c(3,5,2,5))
    if (!is.null(text)) {graphics::par(oma=c(0,0,1,0))}
    
    graphics::plot(TS$Date, TS$Flow, pch=19, cex=0.3, ylab=y1, xlab="",
         xlim=myxlims, ylim=myylims, col="grey50")
    graphics::title(main=mtitle)
    
    if (!is.null(text)) {
        
        if (text == "d") {
            stinfo <- station.info(TS, Plot=F)
            text <- paste("ID: ", stinfo[1],", NAME: ", stinfo[2],
                          ", PROV/STATE: ", stinfo[3], sep = "")
        }
        
        graphics::mtext(text, side=3, line=0, outer=T, cex=0.7)
    }

    ### add shaded polygons covering missing data periods
    MissingDays <- NA.runs(TS)
    polycol <- grDevices::gray(0.5, alpha=0.5)
    if (length(MissingDays$Start) != 0) {
        for (i in 1:length(MissingDays$Start)) {
            xx <- c(MissingDays$Start[i], MissingDays$Start[i],
                    MissingDays$End[i], MissingDays$End[i])
            yy <- c(-100, 1.5*max(TS$Flow), 1.5*max(TS$Flow), -100)
            graphics::polygon(xx,yy,col=polycol, border=NA)
        }
    }
    
    CountsH <- CountsH[CountsH < YearLast - 1]
    CountsB <- CountsB[CountsB < YearLast - 1]
    CountsL <- CountsL[CountsL < YearLast - 1]
    
    hv1 <- graphics::hist(CountsH, breaks=Year_List, plot=F)$counts
    hv2 <- graphics::hist(CountsL, breaks=Year_List, plot=F)$counts
    hv3 <- graphics::hist(CountsB, breaks=Year_List, plot=F)$counts
    
    ## add histogram of changepoints
    if (type == "a") {
        graphics::par(new = T)
        graphics::barplot(rbind(hv3, hv2, hv1), col=mcol, xlab="", ylim=c(0, 30), 
                xaxt="n", yaxt="n")
    }
    if (type == "h") {
        graphics::par(new = T)
        graphics::barplot(hv1, col=mcol, xlab="", ylim=c(0, 10), 
                xaxt="n", yaxt="n")
    }
    if (type == "l") {
        graphics::par(new = T)
        graphics::barplot(hv2, col=mcol, xlab="", ylim=c(0, 10), 
                xaxt="n", yaxt="n")
    }
    if (type == "b") {
        graphics::par(new = T)
        graphics::barplot(hv3, col=mcol, xlab="", ylim=c(0, 10), 
                xaxt="n", yaxt="n")
    }
    
    graphics::axis(side=4)
    graphics::mtext(side = 4, line=3, rylab)
    
    
    graphics::legend("topright", legend = c("Missing Data Periods", types), col=c(polycol, mcol)
           , fill=c(polycol, mcol), bty="n")
    
    Year_List <- Year_List[1:(NumYears-1)]
    res <- data.frame(Year=Year_List, HighFlow=hv1, LowFlow=hv2, Baseflow=hv3)
    
    if (type == "a") {return(res)}
    
    # restore settings on exit
    on.exit(suppressWarnings(graphics::par(opar)))
    
}


