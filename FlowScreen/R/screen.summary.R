#' Create a screening plot
#'
#' Produces summary screening plots of high flow, low flow, or baseflow metrics.  
#' Each plot shows significant temporal trends and step changes. Intended for use as 
#' a data quality screening tool aimed at identifying streamflow records with 
#' anthropogenic impacts or data inhomogeneities.
#' @param metrics output from \code{\link{metrics.all}}
#' @param type Character indicating the set of metrics to plot.  Options
#'   are "h" for high flow metrics, "l" for low flow metrics, or "b" for baseflow
#'   metrics.
#' @param language Language for plot labels.  Choice of either "English" or
#'   "French". Default is "English".
#' @param StnInfo Optional data.frame containing user-supplied station info for plot. 
#'   data.frame must have 7 columns containing station info in the following order:
#'   Station ID, Station Name, Prov/State, Country, Latitude, Longitude, Catchment Area
#'   If any of the information is unavailabe, fill with NA.  The Station ID column must
#'   match the Station ID in column 1 of the data.frame input from \code{\link{create.ts}}.
#' @details For the center of volume (COV) plots on the high flow and baseflow screening plots,
#'   the correlation coefficients for COV and years and for
#'   mean annual flow (MAF) and years are added to the plot. The ratio of the correlation coefficients 
#'   (r COV-years / r COV-MAF) is included as a rudimentary indication
#'   of whether or not the temporal trend in COV is meaningful. See Whitfield (2013) for a
#'   discussion of COV.
#'   
#'   Important note:  If "French" is the language wanted for the plot labels, the language
#'   option must also be specified in \code{\link{metrics.all}}, as this plotting function 
#'   pulls the metric names from the output metrics.all output.
#' @references Whitfield, P.H. 2013. Is 'Center of Volume' a robust indicator of changes 
#'   in snowmelt timing? Hydrological Processes 27:2691-8.
#' @author Jennifer Dierauer
#' @export
#' @examples
#' # load results from metrics.all function for the Caniapiscau River
#' data(caniapiscau.res)
#' 
#' # create a summary flow screening plot of the high flow metrics
#' screen.summary(caniapiscau.res, type="l")


screen.summary <- function(metrics, type="h", language="English", StnInfo=NULL) {
    
    opar <- graphics::par(no.readonly = TRUE)
    
    TS <- metrics[[3]]
    inmetrics <- metrics[[1]]
    inparams <- metrics[[2]]
    
    # calculate mean annual flow correlation with time for cov plots
    maf <- inmetrics[[21]]

    ## set up vector for symbol colors based on data types
    if (type=="l") {DataType <- c(1,1,2,2,2,3,5,1,1,1)} # 5 = km3
    if (type=="h") {DataType <- c(1,2,1,3,1,1,2,2,2,3)} # Q = 1, DOY = 2, #days=3
    if (type=="b") {DataType <- c(1,5,1,1,1,4,2,2,2,3)} #4 % and BFI

    # list positions to subset input metrics time series and trend + changepoint info
    if (type == "h") {start <- 1; end <- 10}
    if (type == "l") {start <- 11; end <- 20}
    if (type == "b") {start <- 21; end <- 30}
    
    inmetrics <- inmetrics[start:end]
    inparams <- inparams[start:end]
    
    # get axis lables in correct language
    Qmax <- 0.95
    if (type == "h") {
        temp <- inparams[[3]]$MetricName
        msplit <- strsplit(temp, " ", fixed=T)[[1]]
        Qmax <- substr(msplit[4], 3, 4)
        Qmax <- paste("0.", Qmax, sep="")
    }
    
    MyXlabs <- get.titles.internal(type, language, Qmax)$Xlabs
    MyYlabs <- get.titles.internal(type, language, Qmax)$Ylabs

    ###layout map sheet
    graphics::layout(matrix(c(1,2,3,3,4,5,3,3,6,7,8,9,10,11,12,13),4,4,byrow=TRUE))
    graphics::par(oma=c(0,0,3,0))
    
    # check formatting of input StnInfo, if StnInfo is not NULL
    if (!is.null(StnInfo)) {
        StnInfo[,1] <- as.character(StnInfo[,1])
        StnInfo[,2] <- as.character(StnInfo[,2])
        StnInfo[,3] <- as.character(StnInfo[,3])
        StnInfo[,4] <- as.character(StnInfo[,4])
        StnInfo[,5] <- as.numeric(StnInfo[,5])
        StnInfo[,6] <- as.numeric(StnInfo[,6])
        StnInfo[,7] <- as.numeric(StnInfo[,7])
    }

    ### Panel 1 - Station Info Text #####
    StnInfo <- station.info(TS, StnInfo, Plot=T, language)
    Country <- StnInfo$Country

    ### Plot 2 - flow regime summary plot
    regime.internal(TS)

    ### Add Main Plot Title
    Agency <- TS$Agency[1]
    if (is.na(Agency)) {Agency <- "Unknown"}
    if (is.na(Country)) {Country <- "Unknown"}
    
    if (language == "English") {
        tstid <- "Station ID: "
        tagency <- ", Agency: "
        tcountry <- ", Country: "
    }
    if (language == "French") {
        tstid <- "ID station: "
        tagency <- ", Agence: "
        tcountry <- ", Pays: "
    }
    graphics::mtext(paste(tstid, TS$ID[1], tagency, Agency, tcountry, Country, sep=""),
              side=3, line=1, outer=TRUE, cex=0.9)


    ### Panel 3 - Larger Q vs time plot
    if (Agency == "WSC") {
        
        SYMs <- c("", "E", "A", "B", "D", "R")
        SYMnames <- c("No Code", "Estimate", "Partial Day", "Ice Conditions", "Dry", "Revised")
        SYMcol <- c("black", "#E41A1C", "#4DAF4A", "#377EB8", "#FF7F00", "#984EA3")
        codes <- as.factor(TS$Code)
        codes <- match(codes, SYMs)
        graphics::par(mar=c(4,4,0,0.5))
        mYlims <- c(0, 1.2*max(TS$Flow))
        graphics::plot(TS$Date, TS$Flow,
             pch=19, col=SYMcol[codes], cex=0.5,
             xlab="", ylab="", ylim=mYlims)
        graphics::title(ylab=MyYlabs[3], line=2)
        graphics::legend(TS$Date[1], 1.15*max(TS$Flow), SYMnames, pch=19, pt.cex=0.9, cex=0.9, col=SYMcol,
               bty="n", xjust=0, x.intersp=0.5, yjust=0.5, ncol=3)
        
    } else {
        
        graphics::par(mar=c(4,4,0,0.5))
        mYlims <- c(0, 1.2*max(TS$Flow))
        graphics::plot(TS$Date, TS$Flow,
             pch=19, cex=0.5,
             xlab="", ylab="", ylim=mYlims)
        graphics::title(ylab=MyYlabs[3], line=2)
        
    }

    ### add shaded polygons covering missing data periods
    MissingDays <- NA.runs(TS)
    polycol <- grDevices::heat.colors(1, alpha=0.2)
    if (length(MissingDays$Start) != 0) {
        for (i in 1:length(MissingDays$Start)) {
            xx <- c(MissingDays$Start[i], MissingDays$Start[i],
                    MissingDays$End[i], MissingDays$End[i])
            yy <- c(-100, 1.5*max(TS$Flow), 1.5*max(TS$Flow), -100)
            graphics::polygon(xx,yy,col=polycol, border=NA)
        }
    }
        
    
    Year1 <- min(c(as.numeric(TS$hyear[1]), as.numeric(TS$year[1])))
    YearEnd <- max(c(max(as.numeric(TS$year)), max(as.numeric(TS$hyear))))
    hyrstart <- as.numeric(subset(TS, TS$hmonth==1)$month[1])
    
    ## Plots 4 to 13
    for (i in 1:10) {
        
        screen.summary.internal(inmetrics[[i]], inparams[[i]], type, 
                                   MyYlabs, i, DataType, maf, Year1, YearEnd, hyrstart)
    }

    on.exit(suppressWarnings(graphics::par(opar)))
    on.exit(graphics::layout(1,1,1))
}
