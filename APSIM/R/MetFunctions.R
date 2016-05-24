#' An S4 class used to store all information about a met file.
#' 
#' @export
#' @slot const A character vector containing constants that wilkl be written to the file. Format is 'variable = value'.
#' @slot lat A length one numeric vector.
#' @slot lon A length one numeric vector.
#' @slot tav A length one numeric vector.
#' @slot amp A length one numeric vector.
#' @slot units A character vector containing the names of the columns in \code{data}.
#' @slot data A data frame containing per day weather data.
metFile <- methods::setClass("metFile",
                    slots = c(const="character", lat="numeric", lon="numeric", tav="numeric", amp="numeric", units="character", data="data.frame"))

#' Convert raw data to the correct APSIM met format.
#' 
#' \code{prepareMet} accepts a data frame containing met data and prepares it 
#' for writing to an APSIM formatted source file.
#' 
#' It will generate year/day columns from an existing date column (and add the required units), check to 
#' ensure that dates are continuous and checks for the existence of required 
#' column names (year, day, radn, mint, maxt and rain).
#' 
#' It will interpolate 365 day leap years (e.g. no extra day from GCMs) and 
#' returns a metFile object that can be used with other APSIM functions.
#' @section Importing External Data: \code{prepareMet} accepts a standard R data
#'   frame as an argument. As such, you can use any importation package that 
#'   returns data in or can be coerced to a data frame. Some examples: #' 
#'   \itemize{ \item Microsoft Excel files - readxl \item NetCDF - RNetCDF \item
#'   MySQL database - RMySQL \item Generic databases (including Microsoft SQL 
#'   Server) - RODBC }
#'   
#' @section Specifiying units for data: Each data column requires a unit in order to be valid. Units need to be enclosed in parentheses.
#'   For unitless values, use "()". Units can be specified by passing a 'units' vector to \code{prepareMet} (see example) or they may already be 
#'   included in the data as would be seen in an APSIM output file. In this case, use \code{\link{loadMet}} instead.
#' @param data A data frame containing the data to prepare.
#' @param lat Latitude in decimal degrees.
#' @param lon Longitude in decimal degrees.
#' @param units A character vector containing units for each column.
#' @param newNames (optional) A vector of new column names.
#' @param date.format (optional) A string containg the date format to use.
#' @return A metFile S4 class containing the prepared met data.
#' @export
#' @examples
#' data(Kingsthorpe)
#' newNames <-c("Date", "maxt", "mint", "rain", "evaporation", "radn", "vp", "Wind", "RH", "SVP")
#' units <- c("()", "(oC)", "(oC)", "(mm)", "(mm)", "(MJ/m^2/day)", "()", "()", "()", "()")
#' prepareMet(kingsData, -27.48, 151.81, newNames = newNames, units = units)
prepareMet <- function (data, lat=stop("Latitude required."), lon=stop("Longitude required."), units=stop("Vector for met units required. See function help if part of table."), newNames=NULL, date.format="AU"){
    #rename columns
    if (!is.null(newNames) & (length(newNames) == length(data))) {
        names(data) <- newNames
    }
    
    # check for existance of year and day
    colNames <- names(data)
    reqNames <- c("year", "day")
    
    if (!all(reqNames %in% colNames)) {
        #at least one year/day column is missing so look for a date column
        
        ifelse(date.format == "AU", date.format <- "%d/%m/%Y",
        ifelse(date.format == "DMY.", date.format <- "%d.%m.%Y",
        ifelse(date.format == "MDY.", date.format <- "%m.%d.%Y",
        ifelse(date.format == "YMD", date.format <- "%Y-%m-%d",
        ifelse(date.format == "US", date.format <- "%m/%d/%Y", date.format <- date.format)))))
        
        converted <- FALSE
        for (col.idx in seq_len(ncol(data))) {
            x <- data[, col.idx]
            if (lubridate::is.Date(x)){
                 converted = TRUE
                 break
            }
            
            if (!is.character(x) | is.factor(x)) next
            if (all(is.na(x))) next
            
                complete.x <- !(is.na(x))
                d <- as.Date(lubridate::parse_date_time(as.character(x), date.format, quiet = TRUE))
                d.na <- d[complete.x]
                if (any(is.na(d.na))) next
                data[, col.idx] <- d
                converted <- TRUE
            
            if (converted) break
        }
        
        if(converted) {
            # found a date column; turn it into year and day
            data$year <- lubridate::year(data[, col.idx])
            data$day  <- lubridate::yday(data[, col.idx])
            units <- c(units, "()", "()")
            colnames(data)[col.idx] <- "date"
        }
        else stop(paste("Could not find year/day or date columns or date column format does not match", date.format))
    }
    
    # check for continuity in dates
    plyr::ddply(data, "year", checkCont)
    
    # Check for standard met columns
    reqNames <- c("maxt", "mint", "radn", "rain", "year", "day")
    print("Required column name check:")
    print(reqNames)
    print(reqNames %in% names(data))
    if (!all(reqNames %in% names(data))) stop("One or more required column names are missing.")
    
    if(ncol(data) != length(units)) stop("The number of columns in the data did not match the number of units. All data columns must have units. For unitless values use ().")
    
    met <- metFile(data=data, lat=lat, lon=lon, units= units)
    
    # add tav and amp
    met <- insertTavAmp(met)
    
    return(met)
}

#' Check for met file errors.
#' 
#' Checks for errors as described in
#' Wall, B.H. "TAMET: Computer program for processing meteorological data." CSIRO
#' Australia. Division of Tropical Crops and Pastures.Tropical Agronomy Technical Memorandum
#' (1977): No. 4, 13p.
#' 
#' Errors checked include:
#' \itemize{
#'   \item Temperature discontinuites.
#'   \item Temperatures that are too high orlow.
#'   \item Evaporation that is too high or low.
#'   \item Radiation that is too high or low.
#'   }
#'   
#'   Note that issues found may not stop APSIM from running
#'   but might indicate an issue with the weather data.
#'   Warnings may not be applicable for very hot or cold climates.
#' 
#' Expects input in metFile object. Use prepareMet or loadMet first.
#'   
#' @param met Met file object.
#' @param lmint Lower bound on minimum temperature.
#' @param umint Upper bound on minimum temperature.
#' @param lmaxt Lower bound on maximum temperature.
#' @param umaxt Upper bound on maximum temperature.
#' @export
#' @examples
#' data(met)
#' checkMet(met)
checkMet <- function (met, lmint=-8, umint=32, lmaxt=10, umaxt=50){
    if (class(met) != "metFile") stop ("Error: checkMet expects a metFile object. Have you run prepareMet?")
    
    if (length(met@lat) == 0 || abs(met@lat) > 90) stop ("latitude not given or out of range.")
    
    met@data$maxtP1 <- c(rep(NA,1),head(met@data$maxt,-1))
    met@data$maxtP2 <- c(rep(NA,2),head(met@data$maxt,-2))
    met@data$mintP1 <- c(rep(NA,1),head(met@data$mint,-1))
    met@data$mintP2 <- c(rep(NA,2),head(met@data$mint,-2))
    
    # get maximum extraterrestrial radiation
    rcal <- sirad::extrat(met@data$day, sirad::radians(met@lat))
    
    # extract evaporation if it exists
    if("evap" %in% names(met@data)){
        evapCol <- met@data[,c("evap", "year", "day")]
    } else if("evaporation" %in% names(met@data)) {
        evapCol <- met@data[,c("evaporation", "year","day")]
        } else {
            evapCol <- NULL
        }
    
    # do some more checks. Not the fastest given the loop, but this isn't run often.
    for(i in 1:nrow(met@data)) {
        # Check for maxt discontinuity.
        if(i > 3 && abs(met@data$maxtP1[i] - met@data$maxt[i] + met@data$maxtP1[i] - met@data$maxtP2[i]) > (ifelse(abs(met@lat) > 18, 19, 9)))
            warning(paste("Maximum temperature discontinuity found around", met@data$year[i], met@data$day[i]))
        
        # Check for mint discontinuity.
        if(i > 3 && abs(met@data$mintP1[i] - met@data$mint[i] + met@data$mintP1[i] - met@data$mintP2[i]) > 20)
            warning(paste("Minimum temperature discontinuity found around", met@data$year[i], met@data$day[i]))
        
        # Evaporation (if available) must be between -0.5 and 20 mm
        # first check for an evaporation column
        if(class(evapCol) != "NULL")
            if (evapCol[i,1] > 20 | evapCol[i,1] < -0.5)
                warning(paste("Evaporation out of bounds on", evapCol$year[i], evapCol$day[i]))
        
        # maxt must be between 10 and 50 degrees Celcius
        if(met@data$maxt[i] < 10 | met@data$maxt[i] > 50)
            warning(paste("Maximum temperature", met@data$maxt[i],"out of range for", met@data$year[i], met@data$day[i]))
        
        # mint must be between -8 and 32 degrees Celcius
        if(met@data$mint[i] < -8 | met@data$mint[i] > 32)
            warning(paste("Minimum temperature", met@data$mint[i]," out of range for", met@data$year[i], met@data$day[i]))
        
        # check for very low radiation values
        if(met@data$radn[i] < rcal$ExtraTerrestrialSolarRadiationDaily[i] * 0.75 * 0.1)
            warning(paste("Radiation value of", met@data$radn[i],"is low for", met@data$year[i], met@data$day[i]))
        
        # check for high radiation values
        if(met@data$radn[i] > rcal$ExtraTerrestrialSolarRadiationDaily[i] * 0.75 * 1.12)
               warning(paste("Radiation value of", met@data$radn[i],"is high for", met@data$year[i], met@data$day[i]))
    }
}

#' Checks a single year for continuity. Called from prepareMet.
#' 
#' This is an internal function used by prepareMet to check for
#' continuity in a single year of a met file.
#' 
#' From tav_amp.for in the APSIM source code:
#'      One earth orbit around the sun does not take an integral
#'      number of days - 365 + a small part of a day.  Since the
#'      Gregorian calendar year is measured as 365 days, a correction
#'      for this err is made every fourth year by adding one day
#'      to the length of the year.  This correction is a little too
#'      much, thus in the centesimal years the correction is not made.
#'      However this over corrects, so in every fourth centesimal year,
#'      the correction of adding one day is made.
#'
#'          To summarise -
#'          If the year is divisible by 4 it is a leap year, unless it is
#'      a centesimal year, in which case it must be divisible by 400.
#'      i.e.  it is a leap year if either of the conditions hold:
#'             (1) the year is divisible by 4 but not by 100;
#'             (2) the year is divisible by 400.
#' @param data A dataframe containing met file data.
#' @export
checkCont <- function(data){
    if((data$year[1] %% 4 == 0 && data$year[1] %% 100 == 0) || data$year[1] %% 400 == 0){
        if(nrow(data) == 365){
            # we're missing a leap day. Interpolate and add extra day
            dr <- data[59,]
            dr$day <- dr$day + 1
            dr$date <- dr$date + 1
            dr$mint <- data$mint[59] + data$mint[60] / 2
            dr$maxt <- data$maxt[59] + data$maxt[60] / 2
            dr$radn <- data$radn[59] + data$radn[60] / 2
            dr$rain <- 0 # rain is too variable to interpolate so just set to 0
            data$day[60:365] <- data$day[60:365] + 1
            data$date[60:365] <- data$date[60:365] + 1
            data <- rbind(data[1:59,], dr, data[60:nrow(data-60),])
        }        
    }   
    
    if(nrow(data) != data$day[nrow(data)] - data$day[1] + 1)
           warning(paste("Number of days in", data$year[1], "does not match expected number,", data$day[nrow(data)] - data$day[1] + 1))
}

#' Inserts Tav and Amp into a met object.
#' 
#' Amp is obtained by averaging the mean daily temperature of each month over 
#' the entire data period resulting in twelve mean temperatures, and then 
#' subtracting the minimum of these values from the maximum. Tav is obtained by 
#' averaging the twelve mean monthly temperatures.
#' 
#' The original documentation for the stand alone Tav_Amp program can be found at
#' \url{http://www.apsim.info/Portals/0/OtherProducts/tav_amp.pdf}.
#' @param met A met file object where the tav and amp will be inserted.
#' @return A met file object to which tav and amp has been added.
#' @export
#' @examples
#' data(met)
#' insertTavAmp(met)
insertTavAmp <- function(met){
    data <- met@data
    #add a month column
    data$month <- lubridate::month(as.Date(paste(data$year, data$day,sep="-"), format="%Y-%j"))
    data$meanDayT <- (data$maxt + data$mint) / 2
    mmt <- plyr::ddply(data, "month", function(df) mean(df$meanDayT))
    if (nrow(mmt) != 12) stop("At least 12 months of data is required to generate tav and amp.")
    met@tav <- max(mmt$V1) - min(mmt$V1)
    met@amp <- mean(mmt$V1)
    return(met)
}

#' Read an APSIM formatted met file.
#' 
#' This function reads in a previously formatted met file. Can be useful for
#' calculating tav and amp for a previously completed file or for reading into R
#' for further analysis.
#' 
#' For reading in other formats see the Importing External Data section in \code{\link{prepareMet}}.
#' @param f The full path to the file to read.
#' @return A metFile object containing the read met file.
#' @export
#' @examples
#' \dontrun{loadMet("Weather.met")}
loadMet <- function(f)
{
    con <- file(f, open="r")
    latFound <- FALSE
    lonFound <- FALSE
    tavFound <- FALSE
    ampFound <- FALSE
    namesFound <- FALSE
    unitsFound <- FALSE
    constants <- NULL
    met <- metFile()
    count <- 0
    data <- list(NULL)
    size <- 1    
    
    while (length(oneLine <- readLines(con, n=1, warn=FALSE)) > 0) {
        #clear out any extra white space
        oneLine <- stringr::str_trim(oneLine)
        
        # skip lines starting with a comment or '[weather' or are blank
        if(ifelse(is.na(stringr::str_locate(oneLine, "!")[1]), FALSE, stringr::str_locate(oneLine, "!")[1] == 1) || 
           ifelse(is.na(stringr::str_locate(oneLine, stringr::fixed("[weather"))[1]), FALSE, stringr::str_locate(oneLine, stringr::fixed("[weather"))[1] == 1) ||
           nchar(oneLine) == 0) next
        
        # look for values
        if(grepl("lat", tolower(oneLine)) && !latFound) {               # look for a latitude
            met@lat <- as.numeric(stringr::str_extract(oneLine, "[-+]?[0-9]*\\.?[0-9]+"))
            latFound <- TRUE
        }else if(grepl("lon", tolower(oneLine)) && !lonFound) {         # look for a longitude
            met@lon <- as.numeric(stringr::str_extract(oneLine, "[-+]?[0-9]*\\.?[0-9]+"))
            lonFound <- TRUE
        }else if(grepl("tav", tolower(oneLine)) && !tavFound) {         # look for a tav
            met@tav <- as.numeric(stringr::str_extract(oneLine, "[-+]?[0-9]*\\.?[0-9]+"))
            tavFound <- TRUE
        } else if(grepl("amp", tolower(oneLine)) && !ampFound) {        # look for an amp
            met@amp <- as.numeric(stringr::str_extract(oneLine, "[-+]?[0-9]*\\.?[0-9]+"))
            ampFound <- TRUE
        } else if(grepl("=", oneLine)) {     #add constants
            constants[length(constants) + 1] <- oneLine
        } else { # we're up to the data now
            if (!namesFound) {
                colNames <- unlist(strsplit(oneLine, " ", fixed=TRUE))
                colNames <- subset(colNames, colNames != "")
                namesFound <- TRUE
            } else if (!unitsFound) {
                units <- unlist(strsplit(oneLine, " ", fixed=TRUE))
                units <- units[units != ""]
                met@units <- units
                if(length(units) != length(colNames)) stop(paste("Error reading", f, "number of columns does match number of headings."))
                unitsFound <- TRUE
            } else {                
                if(count == size) {
                    length(data) <- size <- size * 2
                }
                count <- count + 1
                dataLine <- unlist(strsplit(oneLine, "\\s"))
                dataLine <- dataLine[dataLine != ""]
                data[[count]] <-  dataLine
            }
        }
    }
    close(con)
    data <- data[!sapply(data, is.null)]
    data <- data.frame(matrix(unlist(data), nrow=length(data), byrow=T), stringsAsFactors=FALSE)
    names(data) <- colNames
    #convert character columns to numeric where possible
    for (i in 1:ncol(data)){
        if(!any(is.na(as.numeric(data[[i]]))))
            data[[i]] <- as.numeric(data[[i]])
    }
    met@data <- data
    met@const <- ifelse(is.null(constants), "", constants)
    return(met)
}

#' Write a metFile object to disk in APSIM met format.
#' 
#' Takes a completed metFile object (generated via \code{\link{prepareMet}} or
#' \code{\link{loadMet}}) and writes it to disk.
#' 
#' Note the file will not be written if the met object is missing required
#' information such as latitude, tav or amp. Note that while longitude is required
#' in \code{\link{prepareMet}}, it is not strictly required by APSIM and thus is
#' optional here. However, it is best practise to include it so it will remain
#' mandatory when building a met file via \code{\link{prepareMet}}.
#' @param fileName The file name to write to.
#' @param met The metFile object to write.
#' @export
writeMetFile <- function(fileName, met){
    # do some checks
    if(length(met@lat) == 0) stop("No latitude found.")
    if(length(met@tav) == 0) stop("No tav found. Run insertTavAmp.")
    if(length(met@amp) == 0) stop("No amp found. Run insertTavAmp.")
    
    con <- file(fileName, "w")
    writeLines("[weather.met.weather]", con)
    for(i in 1:length(met@const)) {
        writeLines(met@const[i], con)
    }
    writeLines("", con)
    writeLines(paste("Latitude =", met@lat), con)
    writeLines("", con)
    if(length(met@lon) > 0){
        writeLines(paste("Longitude =", met@lon), con)
        writeLines("", con)
    }
    writeLines(paste("tav =", met@tav), con)
    writeLines(paste("amp =", met@amp), con)
    writeLines("", con)
    writeLines(paste(names(met@data), sep="", collapse=" "), con)
    writeLines(paste(met@units, sep="", collapse=" "), con)
    write.table(met@data, file=con, quote=FALSE, row.names=FALSE, col.names=FALSE)
    close(con)
}