#' Create a plot of the daily streamflow time series
#'
#' Plots the daily streamflow time series and color codes points by data quality codes if 
#' the data are from Water Survey Canada. Also highlights date ranges with missing
#' observations.
#' @param TS output from \code{\link{create.ts}} containing a data.frame of flow
#'   time series
#' @param StnInfo Optional data.frame containing user-supplied station info for plot. 
#'   data.frame must have 7 columns containing station info in the following order:
#'   Station ID, Station Name, Prov/State, Country, Latitude, Longitude, Catchment Area
#'   If any of the information is unavailabe, fill with NA.  The Station ID column must
#'   match the Station ID in column 1 of the data.frame input from \code{\link{create.ts}}.
#' @param text optional character string for margin text, e.g. for station name, 
#'   location, or other notes. Set to NULL if not margin text is wanted, or set to "d" 
#'   to use default text containing the station ID, station name, and province/state 
#'   returned from \code{\link{station.info}}. 
#' @author Jennifer Dierauer and Paul Whitfield
#' @export
#' @examples
#' # load flow time series for the Caniapiscau River
#' data(cania.sub.ts)
#' 
#' # plot daily time series with default margin text
#' screen.series(cania.sub.ts)

screen.series <- function (TS, StnInfo = NULL, text="d") {
    
    opar <- graphics::par(no.readonly = TRUE)

    y1 <- expression(paste("Discharge (m" ^{3}, "/s)"))
    
    if (!is.null(text)) {graphics::par(oma = c(0, 0, 3, 0))}
    
    stdata <- station.info(TS, StnInfo)
    Country <- stdata[4]
    
    Agency <- TS$Agency[1]
    if (is.na(Agency)) {
        Agency <- "Unknown"
    }
    if (is.na(Country)) {
        Country <- "Unknown"
    }
    
    if (Agency == "WSC") {
        SYMs <- c("", "E", "A", "B", "D", "R")
        SYMnames <- c("No Code", "Estimate", "Partial Day", "Ice Conditions", 
                      "Dry", "Revised")
        SYMcol <- c("black", "#E41A1C", "#4DAF4A", "#377EB8", 
                    "#FF7F00", "#984EA3")
        codes <- as.factor(TS$Code)
        codes <- match(codes, SYMs)
        graphics::par(mar = c(4, 3, 0, 0.5))
        mYlims <- c(0, 1.2 * max(TS$Flow))
        graphics::plot(TS$Date, TS$Flow, pch = 19, col = SYMcol[codes], 
             cex = 0.5, xlab = "", ylab = "", ylim = mYlims)
        graphics::title(ylab = y1, line = 2)
        graphics::legend(TS$Date[1], 1.15 * max(TS$Flow), SYMnames, pch = 19, 
               pt.cex = 0.9, cex = 0.9, col = SYMcol, bty = "n", 
               xjust = 0, x.intersp = 0.5, yjust = 0.5, ncol = 3)
    }
    else {
        graphics::par(mar = c(3, 4, 0, 0.5))
        mYlims <- c(0, 1.2 * max(TS$Flow))
        graphics::plot(TS$Date, TS$Flow, pch = 19, cex = 0.5, xlab = "", ylab="", ylim = mYlims)
        graphics::title(ylab = y1, line = 2)
    }
    
    MissingDays <- NA.runs(TS)
    polycol <- grDevices::heat.colors(1, alpha = 0.2)
    if (length(MissingDays$Start) != 0) {
        for (i in 1:length(MissingDays$Start)) {
            xx <- c(MissingDays$Start[i], MissingDays$Start[i], 
                    MissingDays$End[i], MissingDays$End[i])
            yy <- c(-100, 1.5 * max(TS$Flow), 1.5 * max(TS$Flow), 
                    -100)
            graphics::polygon(xx, yy, col = polycol, border = NA)
        }
    }
    
    if (!is.null(text)) {
        
        if (text == "d") {
            stinfo <- station.info(TS, Plot=F)
            text <- paste("ID: ", stinfo[1],", NAME: ", stinfo[2],
                          ", PROV/STATE: ", stinfo[3], sep = "")
        }
        
        graphics::mtext(text, side=3, line=1, outer=T, cex=0.7)
    }
    
    on.exit(suppressWarnings(graphics::par(opar)))
    
}