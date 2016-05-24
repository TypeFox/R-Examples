#' Retrieve Station Info
#'
#' Returns station information for Water Survey Canada or United States Geological Survey stream 
#' gauges.
#' @param TS output from \code{\link{create.ts}} containing a data.frame of flow
#    time series
#' @param StnInfo Optional data.frame containing user-supplied station info for plot. 
#'   data.frame must have 7 columns containing station info in the following order:
#'   Station ID, Station Name, Prov/State, Country, Latitude, Longitude, Catchment Area
#'   If any of the information is unavailabe, fill with NA.  The Station ID column must
#'   match the Station ID in column 1 of the data.frame input from \code{\link{create.ts}}.
#' @param Plot Boolean indicating whether a plot of station information should be created.
#'   Default is F. Plot is intended for use as the upper-left panel of the plot produced 
#'   by \code{\link{screen.summary}}.
#' @param language Language for plotting when Plot = T.  Choice of either "English" or
#'   "French". Default is "English".
#' @return Returns a list of the following station information:
#'   \itemize{
#'     \item $StationID
#'     \item $StnName
#'     \item $Prov/State - Abbreviation for the province or state in which the station is located
#'     \item $Country - Abbreviation for the country in which the station is located
#'     \item $Lat - Latitude of the station
#'     \item $Long - Longitude of the station
#'     \item $Area - Catchment area, in square kilometers
#'     \item $RHN - Boolean indicating whether the station is part of a reference hydrologic network
#'   }
#' @author Jennifer Dierauer
#' @export
#' @examples
#' data(cania.sub.ts)
#' StnInfo <- station.info(cania.sub.ts)

station.info <- function(TS, StnInfo=NULL, Plot=F, language="English") {
    
    if (Plot == T) {
        graphics::par(mar=c(0,0.5,0,0)) # set margins
        graphics::plot(1:10,1:10,pch="",axes=FALSE,ylab="",xlab="") # and plot area
    }
    
    # set language of text
    if (language == "English") {
        tname <- "Name: "
        tprovince <- "Province: "
        tstatus <- "Status: "
        tcatcharea <- "Catchment Area (km2): "
        tperiod <- "Period: "
        tRHBN <- "RHBN Station"
        tRHN <- "RHN Station"
        tProvState <- "Prov/State: "
        tincreasing <- "Blue Trend Line = Increasing Trend"
        tdecreasing <- "Red Trend Line = Decreasing Trend"
        tnoline <- "No Line: p-value > 0.1"
        tthinline <- "0.05 < p-value <= 0.1"
        tmedline <- "0.01 < p-value <= 0.05"
        tthickline <- "p-value <= 0.01"
    }
    
    if (language == "French") {
        tname <- "Nom: "
        tprovince <- "Province: "
        tstatus <- "Statut: "
        tcatcharea <- "Bassin versant (km2): "
        tperiod <- "Periode d'obseration: "
        tRHBN <- "RHBN Station"
        tRHN <- "RHN Station"
        tProvState <- "Province/Etat: "
        tincreasing <- "ligne bleue = tendance a la hausse"
        tdecreasing <- "ligne rouge = tendance a la baisse"
        tnoline <- "pas de ligne: p > 0.1"
        tthinline <- "0.05 < p <= 0.1"
        tmedline <- "0.01 < p <= 0.05"
        tthickline <- "p <= 0.01"
    }
    
    
    Station <- as.character(TS$ID[1]) # get station ID
    Agency <- as.character(TS$Agency[1]) # get agency
    
    if (is.na(Agency)) {Agency <- "Unknown"}
    
    ## check for user-supplied stninfo
    if (is.null(StnInfo)) {
    
        #if no user-supplied station info, and Agency is WSC, get station info from sysdata
        if (Agency =="WSC") {
            
            Year1 <- min(as.numeric(TS$year))
            YearEnd <- max(as.numeric(TS$year))
            
            StationInfo <- get.station.internal(Station)
            
            # compile station info for plot and return
            StnName <- as.character(StationInfo[,2])
            ProvState <- as.character(StationInfo[,4])
            Country <- "CA"
            Lat <- round(as.numeric(as.character(StationInfo[,5])), 2)
            Long <- round(as.numeric(as.character(StationInfo[,6])), 2)
            Area <- StationInfo[,7]

            lname <- StationInfo[18]
            endchar <- substr(lname, nchar(lname), nchar(lname))
            if (identical(endchar, "*")) {
                RHN <- T
            } else (RHN <- F)
            
            out <- list(Station, StnName, ProvState, Country, Lat, Long, Area, RHN)
            
            ### Add Relevant Station info to Plot if Plot == T
            if (Plot == T) {
                graphics::text(1,9.5,paste(tname, StnName, sep=""), adj=c(0,0))
                graphics::text(1,9,paste(tprovince, ProvState,
                               " ", tstatus, as.character(StationInfo[,3]), sep=""), adj=c(0,0))
                graphics::text(1,8.5,paste("Latitude: ", Lat, " ", "Longitude: ", Long, sep=""), adj=c(0,0))
                graphics::text(1,8, paste(tcatcharea, Area, sep=""), adj=c(0,0))
                graphics::text(1,7.5, paste(tperiod, Year1, "-", YearEnd, sep=""),adj=c(0,0))
                
                if (RHN == T) {
                    graphics::text(1, 7, tRHBN, adj=c(0,0))
                }
            }
            
        #if no user-supplied stninfo and agency is not WSC, check for station in USGS list
        } else {
            
            Year1 <- min(as.numeric(TS$year))
            YearEnd <- max(as.numeric(TS$year))
            
            if (Station %in% USGS.site.info$STAID) {
                
                StationInfo <- USGS.site.info[USGS.site.info$STAID==Station,] # get info from sysdata
                StnName <- StationInfo[2][,1]
                Country <- "United States"
                ProvState <- StationInfo[9][,1]
                Lat <- round(as.numeric(StationInfo[7][,1]), 2)
                Long <- round(as.numeric(StationInfo[8][,1]), 2)
                Area <- StationInfo[5][,1]
                if (identical(StationInfo[10][,1], "yes")) {RHN <- T} else {RHN <- F}
                
                out <- list(Station, StnName, ProvState, Country, Lat, Long, Area, RHN)
                
                if (Plot == T) {
                    graphics::text(1, 9.5, paste("Name:", StnName), adj=c(0,0))
                    graphics::text(1, 9, paste(tProvState, ProvState,"   ", "Status: Active", sep=""), 
                         adj=c(0,0))
                    graphics::text(1, 8.5, paste("Latitude: ", Lat, "   ", "Longitude: ", Long, sep=""), adj=c(0,0))
                    graphics::text(1, 8, paste(tcatcharea, Area, sep=""), adj=c(0,0))
                    graphics::text(1, 7.5, paste(tperiod, Year1, "-", YearEnd, sep=""), adj=c(0,0))
                    graphics::text(1, 7, tRHN, adj=c(0,0))
                }
            
            ## if not a WSC or USGS station and no user supplied station info - fill with NA
            } else {
                
                if (Plot == T) {
                    graphics::text(1, 9.5, "Name: NA", adj=c(0,0))
                    graphics::text(1, 9, "Prov/State: NA", adj=c(0,0))
                    graphics::text(1, 8.5, paste("Latitude: NA", "   ", "Longitude: NA"), adj=c(0,0))
                    graphics::text(1, 8, "Catchment Area (km2): NA", adj=c(0,0))
                    graphics::text(1, 7.5, paste(tperiod, Year1, "-", YearEnd, sep=""), adj=c(0,0))
                }
                
                Country <- "Unknown"
                StnName <- "Unknown"
                RHN <- "Unknown"
                
                out <- list(Station, NA, NA, NA, NA, NA, NA, RHN)
            }
        }
    ## if user-supplied station info, use for plot
    } else {
        
        Year1 <- min(as.numeric(TS$year))
        YearEnd <- max(as.numeric(TS$year))
        StationInfo <- StnInfo[StnInfo[,1] == Station,] # get user-supplied station info
        StnName <- as.character(StationInfo[2][,1])
        ProvState <- as.character(StationInfo[3][,1])
        Country <- as.character(StationInfo[4][,1])
        Lat <- round(StationInfo[5][,1], 2)
        Long <- round(StationInfo[6][,1], 2)
        Area <- StationInfo[7][,1]
        RHN <- "Unknown"
        
        if (Plot == T) {
            graphics::text(1, 9.5, paste(tname, StnName, sept=""), adj=c(0,0))
            graphics::text(1, 9, paste(tProvState, ProvState, sep=""), adj=c(0,0))
            graphics::text(1, 8.5, paste("Latitude: ", Lat, "   ", "Longitude: ", Long, sep=""), adj=c(0,0))
            graphics::text(1, 8, paste(tcatcharea, Area, sep=""), adj=c(0,0))
            graphics::text(1, 7.5, paste(tperiod, Year1, "-", YearEnd, sep=""), adj=c(0,0))
        }
        
        out <- list(Station, StnName, ProvState, Country, Lat, Long, Area, RHN)
    }
        
    if (Plot == T) {
        ### Add Legend for plot trend lines
        graphics::text(1,6, tincreasing, col="darkblue", adj=c(0,0))
        graphics::text(1,5.5, tdecreasing, col="darkred", adj=c(0,0))
        graphics::text(3,4.5, tnoline, adj=c(0,0.5))
        graphics::text(3,4, tthinline, adj=c(0,0.5))
        graphics::text(3,3.5,tmedline, adj=c(0,0.5))
        graphics::text(3,3, tthickline, adj=c(0,0.5))
        graphics::points(c(1,2), c(4, 4), type="l", lwd=1)
        graphics::points(c(1,2), c(3.5, 3.5), type="l", lwd=2)
        graphics::points(c(1,2), c(3, 3), type="l", lwd=3)
    }

    names(out) <- c("StationID", "StnName", "Prov/State", "Country", "Lat", "Long",
                    "Area", "RHN")
    return(out)
    
}