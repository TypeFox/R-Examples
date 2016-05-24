#' Convert a data frame of non-equal interval continuous observations into
#' equal interval continuous observations
#'
#' \code{SurvExpand} convert a data frame of non-equal interval continuous
#' observations into equal interval continuous observations. This is useful for
#' creating time-interactions with \code{\link{tvc}}.
#' @param data a data frame.
#' @param GroupVar a character string naming the unit grouping variable.
#' @param Time a character string naming the variable with the interval start
#' time.
#' @param Time2 a character string naming the variable with the interval end
#' time.
#' @param event a character string naming the event variable. Note: must be
#' numeric with 0 indicating no event.
#' @param PartialData logical indicating whether or not to keep only the
#' expanded data required to find the Cox partial likelihood.
#' @param messages logical indicating if you want messages returned while the
#' function is working.
#'
#' @return Returns a data frame where observations have been expanded into
#' equally spaced time intervals.
#'
#' @details The function primarily prepares data from the creation of accurate
#' time-interactions with the \code{\link{tvc}} command.
#' Note: the function will work best if your original time intervals are
#' recorded in whole numbers. It also currently does not support repeated
#' events data.
#'
#' @examples
#' # Load Golub & Steunenberg (2007) Data
#' data("GolubEUPData")
#'
#' # Subset PURELY TO SPEED UP THE EXAMPLE
#' GolubEUPData <- GolubEUPData[1:500, ]
#'
#' # Expand data
#' GolubEUPDataExpand <- SurvExpand(GolubEUPData, GroupVar = 'caseno',
#'                        Time = 'begin', Time2 = 'end', event = 'event')
#'
#' @references Gandrud, Christopher. 2015. simPH: An R Package for Illustrating
#' Estimates from Cox Proportional Hazard Models Including for Interactive and
#' Nonlinear Effects. Journal of Statistical Software. 65(3)1-20.
#'
#' @seealso \code{\link{tvc}}
#' @import data.table
#' @importFrom dplyr group_by mutate select distinct
#' @keywords utilities
#' @export

SurvExpand <- function(data, GroupVar, Time, Time2, event, PartialData = TRUE,
                        messages = TRUE){
    # For CRAN
    FT <- UG <- FinalTime <- StartT <- ST <- FinishT <- NULL

    # Reminder
    if (isTRUE(messages)){
        message('\nReminder (1): Currently SurvExpand does not support repeated events data.\n\n')
        message('Reminder (2): Times should be in 1 unit intervals, otherwise this could take awhile.\n\n')
    }
    # Warnings
    if (class(data[, Time]) != 'numeric'){
        if (isTRUE(messages)){
            message(paste0('Converting ', deparse(substitute(Time)),
                    ' to numeric. Things might get wacky. Please check.'))
        }
        data[, Time] <- as.character(data[, Time])
        data[, Time] <- as.numeric(data[, Time])
    }
    if (class(data[, Time2]) != 'numeric'){
        if (isTRUE(messages)){
            message(paste0('Converting ', deparse(substitute(Time2)),
                    ' to numeric. Things might get wacky. Please check.'))
        }
        data[, Time2] <- as.character(data[, Time2])
        data[, Time2] <- as.numeric(data[, Time2])
    }
    if (class(data[, event]) != 'numeric'){
        if (isTRUE(messages)){
            message(paste0('Converting ', deparse(substitute(event)),
                    ' to numeric. Things might get wacky. Please check.'))
        }
        data[, event] <- as.character(data[, event])
        data[, event] <- as.numeric(data[, event])
    }

    # Create full unit-time data set
    Min <- min(data[, Time])
    Max <- max(data[, Time2])
    FullTimes = Min:Max
    FullTimes <- data.frame(FakeID = 1, FT = FullTimes)
    FullTimes <- data.table(FullTimes, key = "FakeID", allow.cartesian = TRUE)

    UniGroup <- unique(data[, GroupVar])
    UniGroup <- data.frame(FakeID = 1, UG = UniGroup)
    UniGroup <- data.table(UniGroup, key = "FakeID", allow.cartesian = TRUE)

    Full <- UniGroup[FullTimes, allow.cartesian = TRUE]
    Full <- Full[, c('UG', 'FT'), with = FALSE]
    Full <- setkey(Full, key = 'UG')

    if (isTRUE(messages)) message('\nExpanding data.\n')

    # Keep only times needed for survival modelling
    End <- sort(unique(data[, Time2]))
    if (isTRUE(PartialData)){
        if (isTRUE(messages)) message('Keeping only needed observations.\n')
        FullSub <- Full[FT %in% End]
        FullSub <- FullSub[order(UG, FT)]
    }
    else {
        FullSub <- Full
    }

    # Find last observation and subset data so that unit-observations greater
    # than this are excluded
    DGroup <- eval(parse(text = paste0('group_by(data, ', GroupVar, ')')))
    DGroup <- eval(parse(text = paste0('dplyr::mutate(DGroup, FinalTime = max(',
                                    Time2, '))')))
    names(DGroup)[names(DGroup) == GroupVar] <- "UG"
    DGroup <- group_by(DGroup, UG)
    DGroup <- select(DGroup, UG, FinalTime)
    DGroup <- distinct(DGroup, FinalTime)
    DGroup <- data.frame(DGroup)

    DGroup <- data.table(DGroup, key = 'UG', allow.cartesian = TRUE)
    FullLast <- FullSub[DGroup, allow.cartesian = TRUE]
    FullLast <- FullLast[FT <= FinalTime]
    FullLast <- FullLast[, c('UG', 'FT', 'FinalTime'), with = FALSE]

    # Create new time interval start variable
    FullLast <- FullLast[, 'ST' := FT - 1]
    FullLast <- setkey(FullLast, key = 'UG')

    # Merge with original data and create new time intervals
    DataMerge <- data
    Names <- c(GroupVar, Time, Time2, event)
    Names2 <- c('UG', 'StartT', 'FinishT', 'Event')
    DataMerge <- MoveFront(DataMerge, Names)
    DataMerge <- setnames(DataMerge, old = 1:4, new = Names2)
    DataMerge <- data.table(DataMerge, key = 'UG', allow.cartisian = TRUE)

    FullComb <- FullLast[DataMerge, allow.cartesian = TRUE]
    FullComb <- FullComb[, !c('allow.cartisian'),
                        with = FALSE]

    if (isTRUE(messages)) message('Doing a final clean up.')

    Out <- FullComb[StartT <= ST & FinishT >= FT ]
    Out <- setkey(Out, NULL)
    Out <- unique(Out)

    # Final clean up
    Out <- data.frame(Out)
    Out$Event[Out$FinalTime != Out$FT] <- 0
    FinalDrop <- c('StartT', 'FinishT', 'FinalTime')
    Out <- Out[, !(names(Out) %in% FinalDrop)]
    Out <- MoveFront(Out, c('UG', 'ST', 'FT', 'Event'))
    colnames(Out)[1:4] <- Names

    return(Out)
}
