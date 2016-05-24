#' Merge a location data set with a dry time (or other stopping) covariate
#' 
#' 
#' The function merges a location data set with a stopping variable data set.
#' 
#' 
#' Simply merges the data frames and interpolates based on the chosen method.
#' Both data frames have to use the same name for the time variable. Also
#' contains \code{stopType} which = "o" if observed or "p" for interpolated.
#' 
#' The merged data is truncated to the first and last time in the location data
#' set. Missing values in the stopping variable data set can be interpolated by
#' replacing them with zeros (full movement) or first replacing with zeros then
#' using a moving average to smooth the data. Only the missing values are then
#' replace with this smoothed data. This allows a smooth transition to full
#' movement.
#' 
#' @param data Location data.
#' @param stopData stopping variable data set.
#' @param Time.name character naming time index variable in both data sets
#' @param interp method of interpolation.
#' @param win window for "ma0" interpolation method.
#' @param constCol columns in \code{data} for which the user would like to be
#' constant, such as id or sex.
#' @return
#' 
#' Merged data.frame with new column from \code{stopData}. Missing values in
#' the stopping variable will be interpolated
#' @author Devin S. Johnson
#' @examples
#' 
#' 
#' track <- data.frame(TimeVar=sort(runif(20,0,20)), x=1:20, y=20:1)
#' track
#' stopData <- data.frame(TimeVar=0:29, stopVar=round(runif(30)))
#' stopData
#' mergeTrackStop(track, stopData, Time.name="TimeVar")
#' 
#' @export
"mergeTrackStop" <- function(data, stopData, Time.name="Time",
                             interp=c('zeros','ma0'), win=2, constCol)
{
    interp <- interp[1]
    nmsStop <- names(stopData)[!names(stopData) %in% Time.name]
    if (length(nmsStop) > 1) stop("You can only merge 1 stopping variable at a time")
    stopData <- stopData[order(stopData[, Time.name], stopData[, nmsStop]), ]
    stopData <- stopData[!duplicated(stopData[, Time.name]), ]
    stopStart <- min(stopData[, Time.name][!is.na(stopData[, nmsStop])])
    stopEnd <- max(stopData[, Time.name])
    trackStart <- min(data[, Time.name])
    trackEnd <- max(data[, Time.name])
    Start <- max(stopStart, trackStart)
    End <- min(stopEnd, trackEnd)
    if (!missing(constCol)) constVal <- data[1, constCol]
    mergeData <- data.frame(seq(floor(Start), ceiling(End), 1))
    names(mergeData) <- Time.name
    stopData <- merge(stopData, mergeData, by=Time.name, all=TRUE)
    stopData$stopType <- ifelse(is.na(stopData[, nmsStop]), "p", "o")
    dtTmp <- ifelse(is.na(stopData[, nmsStop]), 0, stopData[, nmsStop])
    if (interp == "ma0") {
        stopTimeAvg <- filter(dtTmp, rep(1 / (2 * win + 1), 2 * win + 1),
                              "convolution", sides=2)
        dtTmp <- ifelse(is.na(stopData[, nmsStop]), stopTimeAvg, stopData[, nmsStop])
    }
    stopData[, paste(nmsStop, "Orig", sep="")] <- stopData[, nmsStop]
    stopData[, nmsStop] <- ifelse(is.na(dtTmp), 0, dtTmp)
    stopData$locType <- "p"
    data$locType <- "o"
    data <- merge(data, stopData, all=TRUE)
    data <- data[order(data[, Time.name]), ]
    data <- data[!duplicated(data[, Time.name]), ]
    data[, nmsStop] <- round(approx(data[, Time.name], data[, nmsStop],
                                    xout=data[, Time.name], method="constant")$y, digits=3)
    pind <- ifelse(data$stopType == "p", 1, 0)
    pind <- approx(pind, xout=1:length(pind), method="constant", rule=2)$y
    data$stopType <- ifelse(pind == 1, "p", "o")
    if (!missing(constCol)) {
        for (l in 1:length(constCol)) {
            data[, constCol[l]] <- constVal[l]
        }
    }
    out <- data[(data[, Time.name] >= trackStart & data[, Time.name] <= trackEnd), ]
    return(out)
}
