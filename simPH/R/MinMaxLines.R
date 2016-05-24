#' Transform the simulation object to include only the min and max of the
#' constricted intervals, as well as the lower and upper bounds of the middle 50
#' percent of the constricted intervals
#'
#' \code{MinMaxLines} is an internal function to transform the simulation
#' object to include only the min and max of the intervals set by \code{ci} in
#' the \code{coxsim} command, as well as the lower and upper bounds of the
#' middle 50 percent of these intervals. It also returns the medians.
#'
#' @param df a data frame or a simulation class object.
#' @param byVars character vector of the variables to subset the data frame by.
#' The default is \code{'Xj'}.
#' @param hr logical indicating whether or not \code{df} contains a hazard rate.
#' @param strata logical indicating whether or not \code{df} contains a
#' stratified hazard rate.
#' @param clean logical, whether or not to clean up the output data frame to
#' only include \code{byVars}, \code{Min_CI}, \code{Lower50_CI}, \code{median},
#' \code{Upper50_CI}, \code{Max_CI}.
#'
#' @examples
#' # Load Carpenter (2002) data
#' data("CarpenterFdaData")
#'
#' # Load survival package
#' library(survival)
#'
#' # Run basic model
#' M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal +
#'            deathrt1 + acutediz + hosp01  + hhosleng +
#'            mandiz01 + femdiz01 + peddiz01 + orphdum +
#'            vandavg3 + wpnoavg3 + condavg3 + orderent +
#'            stafcder, data = CarpenterFdaData)
#'
#'  # Simulate Hazard Ratios
#'  Sim1 <- coxsimLinear(M1, b = "stafcder",
#'                       Xj = c(1237, 1600),
#'                       Xl = c(1000, 1000),
#'                       qi = "Hazard Ratio",
#'                       spin = TRUE, ci = 0.99)
#'
#' # Find summary statistics of the constricted interval
#' Sum <- MinMaxLines(Sim1, clean = TRUE)
#'
#' @importFrom dplyr group_by_ ungroup mutate distinct_
#' @import lazyeval
#' @importFrom stats median quantile
#' @keywords internals
#' @export

MinMaxLines <- function(df, byVars = "Xj", hr = FALSE, strata = FALSE,
                        clean = FALSE){
    Xj <- QI <- Time <- HRValue <- HRate <- Strata <- NULL

    df <- as.data.frame(df)

    if (isTRUE(hr) & !isTRUE(strata)){
        byVars <- c("Time", "HRValue")
    }
    else if (isTRUE(hr) & !isTRUE(strata)){
        byVars <- c("Time", "HRValue", "Strata")
    }

    df <- group_by_(df, .dots = byVars)

    if (!isTRUE(hr)){
        Linesdf <- mutate(df, Median = median(QI),
                        Max = max(QI),
                        Min = min(QI),
                        Lower50 = quantile(QI, 0.25),
                        Upper50 = quantile(QI, 0.75))

        Linesdf <- distinct_(Linesdf, .dots = byVars, .keep_all = TRUE)
    }
    else if (isTRUE(hr) & !isTRUE(strata)){
        Linesdf <- mutate(df, Median = median(HRate),
                          Max = max(QI),
                          Min = min(QI),
                          Lower50 = quantile(QI, 0.25),
                          Upper50 = quantile(QI, 0.75))

        Linesdf <- distinct_(Linesdf, .dots = c(1, 3), .keep_all = TRUE)
        #Linesdf <- Linesdf[!duplicated(Linesdf[, c(1, 3)]), ]
    }
    else if (isTRUE(hr) & isTRUE(strata)){
        Linesdf <- mutate(df, Median = median(HRate),
                        Max = max(QI),
                        Min = min(QI),
                        Lower50 = quantile(QI, 0.25),
                        Upper50 = quantile(QI, 0.75))

        Linesdf <- distinct(Linesdf, Time, HRValue, Strata, .keep_all = TRUE)
        #Linesdf <- Linesdf[!duplicated(
        #                    Linesdf[, c("Time", "HRValue", "Strata")]), ]
    }

    Linesdf <- ungroup(Linesdf)

    if (isTRUE(clean)){
        Linesdf <- Linesdf[, c(byVars, 'Min', 'Lower50', 'Median', 'Upper50',
                                'Max')]
        names(Linesdf) <- c(byVars, 'Min_CI', 'Lower50_CI','Median',
                            'Upper50_CI', 'Max_CI')
    }
    class(Linesdf) <- 'data.frame'
    return(Linesdf)
}


#' Create a variable of each simulation value's percentile in the distribution.
#'
#' @param SimIn data frame of simulations.
#' @param xaxis character string. The column that will form the x-axis in
#' the plot.
#' @param yaxis character string. The column that will form the y-axis in
#' the plot.
#'
#' @importFrom dplyr group_by mutate
#'
#' @keywords internals
#' @noRd

PercRank <- function(SimIn, xaxis, yaxis = 'QI'){
    Xaxis <- QI <- NULL

    names(SimIn)[names(SimIn) == xaxis] <- 'Xaxis'
    names(SimIn)[names(SimIn) == yaxis] <- 'QI'

    PlainPercRank <- function(x) trunc(rank(x))/length(x)

    Temp <- dplyr::group_by(SimIn, Xaxis)
    Temp <- dplyr::mutate(Temp, PercRank = PlainPercRank(QI))

    # Center on the 50th percentile
    Temp$PercRank <- abs(Temp$PercRank - 0.5)
    Temp$PercRank <- round(abs(0.5 - Temp$PercRank), 1)

    names(Temp)[names(Temp) == 'Xaxis'] <- xaxis
    names(Temp)[names(Temp) == 'QI'] <- yaxis
    Temp
}
