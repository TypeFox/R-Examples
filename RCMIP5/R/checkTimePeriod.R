#' Check for continuous time periods in CMIP5 files
#'
#' Check that all time periods are continuous and present for multi-file ensembles.
#' Before starting to process what may be hundreds or thousands of CMIP5 files,
#' it's a good idea to verify that your file set is complete and not missing any
#' years.
#'
#' @param fileInfo_df data.frame from getFileInfo
#' @return A data frame showing which ensembles are continuous, and which are not.
#' In addition to standard identifying fields in the data frame (domain, model,
#' experiment, variable, and ensemble), this includes:
#'  \item{yrStr}{A concatenation of time strings for all ensembles}
#'  \item{allHere}{a quick check for yr and mon frequency}
#'  \item{startDate}{Earliest (decimal) date for the ensemble}
#'  \item{endDate}{Latest (decimal) date for the ensemble}
#'  \item{file}{The number of files in the ensemble}
#' @details  This function calls \code{\link{getFileInfo}} to scan a directory tree,
#' and then examines the time data in these filenames. These time signatures
#' will be concatenated and an 'allHere' flag returned.
#' @note This only works for files that are in domains 'fx', 'mon', or 'yr'.
#' Decimal time is (year + (month-1)/12).
#' @note Unfortunately it's impossible to automatically check the time signature for
#' sub-monthly frequencies quickly without opening the NetCDF file.
#' @examples
#' \dontrun{
#' checkTimePeriod(getFileInfo())
#' }
#' @seealso \code{\link{getFileInfo}}
#' @export
checkTimePeriod <- function(fileInfo_df) {
    
    # Sanity checks
    stopifnot(is.data.frame(fileInfo_df))
    ddplyFields <- c("domain", "experiment","model","variable","ensemble")
    stopifnot(all(ddplyFields %in% colnames(fileInfo_df)))
    stopifnot("time" %in% colnames(fileInfo_df))
    
    # Use ddply to break up data frame, process and check time field, and return result
    result <- data.frame()
    splitter <- apply( fileInfo_df[, ddplyFields], 1, paste, collapse="-" )
    for(i in unique(splitter)) {
        x <- fileInfo_df[splitter==i,]
        
        # split the time signiture
        curCombo <- matrix(unlist(strsplit(as.character(x$time), '-')),
                           ncol=2, byrow=TRUE)
        
        # Find the starting and ending decimal year
        startYear <- as.numeric(substr(curCombo[,1], 1, 4))
        endYear <- as.numeric(substr(curCombo[,2], 1, 4))
        
        # pull the time step from the domain name
        if(all(x$domain %in% 'fx')) { # fixed
            # don't even process fixed files
            next
        } else if(all(grepl('mon', x$domain))) { # monthly
            timeStep <- 1/12
            # convert to decimal years
            startYear <- startYear+(as.numeric(substr(curCombo[,1], 5, 6))-1)/12
            endYear <- endYear + (as.numeric(substr(curCombo[,2], 5, 6))-1)/12
        } else if(all(grepl('yr', x$domain))) { # annual
            timeStep <- 1
        } else {
            # we can not process sub-monthly time scales because of not standard
            #... year lengths. The user must check the strings by hand.
            timeStep <- NA
        }
        
        if(is.na(timeStep)) {
            allHere <- NA
        } else {
            # Figure out the target date for the start of the next file
            nextYear <- endYear + timeStep
        }
        
        # One file is always complete
        if(length(startYear) == 1) {
            allHere <- TRUE
            # If multiple files, shift indexes to compare the start/stop values
        } else if( !is.na(timeStep) & length(startYear) > 1) {
            startIndex <- c(2:length(startYear))
            endIndex <- c((2:length(startYear))-1)
            allHere <- all(abs(nextYear[endIndex]-startYear[startIndex]) < 1e-6)
        }
        
        # return answering data frame which contains
        #   yrStr - All orginal year strings for reference (useful if something is wrong).
        #   allHere - boolean saying if the strings match up
        #   startDate - earliest time stamp
        #   endDate - latest time stamp
        result <- rbind(result, cbind(
            x[1, ddplyFields],
            data.frame(yrStr=paste(x$time, collapse=';'),
                       allHere=allHere,
                       startDate=min(startYear),
                       endDate=max(endYear),
                       files=length(startYear))
        ))
    }
    result
} # checkTimePeriod
