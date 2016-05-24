#' Calculate daily statistics for dendrometer and environmental data
#'
#' @description The function calculates various daily statistics for dendrometer and environmental data. It either returns multiple statistics for individual sensors, or a single statistic for multiple sensors.
#'
#' @usage daily_stats(dm.data, sensor = 1, value = c("mean", "min",
#'             "max", "sum"), smooth.param = 1)
#'
#' @param dm.data a \code{data.frame} with a timestamp (\code{\%Y-\%m-\%d \%H:\%M:\%S} format) as row names, and dendrometer series in columns. Output as created using code from the \code{Import dendrometer data} vignette, or gap-filled dendrometer series as produced by \code{\link{fill_gaps}}. Environmental data can be specified as well, and should be formatted as dendrometer data.
#' @param sensor a \code{numeric} or \code{character} string specifying the sensor(s) to be used in the function. Defaults to 1 (first column of \code{data.frame}). If "ALL" is specified, a single \code{value} will be calculated or extracted for all series in the \code{data.frame}.
#' @param value a \code{character} string of \code{"mean"}, \code{"min"}, \code{"max"} or \code{"sum"}, specifying the daily statistic to be calculated or extracted. Optional argument for \code{sensor = "ALL"}, defaults to \code{"mean"}. Argument matching is performed.
#' @param smooth.param a \code{numeric} specifying the degree of smoothing. Defaults to 1 (no smoothing). In case smoothing is applied, series should be gap-free or gap-filled.
#'
#' @details The function calculates various daily statistics for dendrometer and environmental data. For \code{sensor} is \code{numeric}, the function returns multiple statistics for a single sensor. For \code{sensor = "ALL"}, the function returns a single statistic (i.e. \code{"mean"}, \code{"min"}, \code{"max"} or \code{"sum"}) for all columns of the \code{data.frame}, whereby \code{"sum"} is particularly relevant for environmental parameters like precipitation.
#'
#' The function includes a smoothing option (argument \code{smooth.param}) particularly for noisy datasets in which outliers may under- or overestimate minimum and maximum stem sizes within days. By default, no smoothing is performed. Smoothing requires gap-free series.
#'
#' @return The function returns:
#'
#' \itemize{\item{for \code{sensor} is \code{numeric}, a \code{data.frame} containing the following columns:}}
#' \item{dmID}{dendrometer ID.}
#' \item{date}{timestamp in \code{\%Y-\%m-\%d} format.}
#' \item{DOY}{day of year.}
#' \item{min}{minimum daily stem size.}
#' \item{mean}{mean daily stem size.}
#' \item{max}{maximum daily stem size.}
#' \item{amplitude}{amplitude of daily stem-size changes (i.e. max - min).}
#' \item{time_min}{timestamp indicating the timing of the minimum.}
#' \item{time_max}{timestamp indicating the timing of the maximum.}
#'
#' \itemize{\item{for \code{sensor} is \code{"ALL"}:}}
#'
#' a \code{data.frame} with a timestamp (\code{\%Y-\%m-\%d}) as row names, and processed dendrometer or environmental data in columns (i.e. mean, minimum, maximum or sum).
#'
#' @author Olivier Bouriaud, Ernst van der Maaten and Marieke van der Maaten-Theunissen.
#'
#' @examples
#' data(dmCD)
#' dm.daily <- daily_stats(dmCD, sensor = 1)
#'
#' data(dmED)
#' dm.daily <- daily_stats(dmED, sensor = "ALL", value = "max")
#'
#' @import stats
#'
#' @export daily_stats
#'
daily_stats <- function(dm.data, sensor = 1, value = c("mean", "min", "max", "sum"), smooth.param = 1)
{
  nm <- deparse(substitute(dm.data))

  if(!is.dendro(dm.data)) {
    stop(paste("'", nm, "' is not in the required format", sep = ""))
  }

  if(class(sensor) == "numeric") {
    if(ncol(dm.data) == 1 && sensor > ncol(dm.data)) {
      stop("'sensor' should be 1")
    }
    if(ncol(dm.data) > 1 && sensor > ncol(dm.data)) {
      stop(paste("'sensor' should be between", 1, "and", ncol(dm.data), sep = " "))
    }

    # Optional smoothing
    if(smooth.param > 1) {
      dm.data[, sensor] <- ave(dm.data[, sensor],
                              FUN = function(x) rollmean(x, smooth.param, align = "right", fill = NA, na.rm = TRUE))
    }

    date <- strftime(rownames(dm.data), format = "%Y-%m-%d")
    DOY <- as.numeric(strftime(rownames(dm.data), format = "%j"))
    dmin <- aggregate(dm.data[,sensor], list(DOY, DOY), FUN = min)$x
    dmean <- aggregate(dm.data[,sensor], list(DOY, DOY), FUN = mean)$x
    dmax <- aggregate(dm.data[,sensor], list(DOY, DOY), FUN = max)$x
    amplitude <- dmax - dmin

    dm.data$DOY <- DOY
    sort.min <- dm.data[order(dm.data$DOY, dm.data[,sensor]),]
    sort.min <- sort.min[!duplicated(sort.min$DOY),]
    time_min <- row.names(sort.min)
    sort.max <- dm.data[order(dm.data$DOY, -dm.data[,sensor]),]
    sort.max <- sort.max[!duplicated(sort.max$DOY),]
    time_max <- row.names(sort.max)

    dmDaily <- data.frame(dmID = rep(names(dm.data)[sensor], length(unique(DOY))),
                           date = unique(date),
                           DOY = unique(DOY),
                           min = dmin,
                           mean = dmean,
                           max = dmax,
                           amplitude = amplitude,
                           time_min = time_min,
                           time_max = time_max)

    dmDaily$time_min[is.na(dmDaily$min)] <- NA
    dmDaily$time_max[is.na(dmDaily$min)] <- NA
  }

  if(class(sensor) != "numeric") {
    if(sensor != "ALL") {
      stop(paste("the entry for 'sensor' is incorrect. Either specify a sensor number or type ", dQuote("ALL"),".", sep = ""))
    }

    date <- strftime(rownames(dm.data), format = "%Y-%m-%d")
    DOY <- as.numeric(strftime(rownames(dm.data), format = "%j"))

    value <- match.arg(value, c("mean", "min", "max", "sum"))

    dmvalue <- aggregate(dm.data, list(DOY, DOY), FUN = value)

    dmDaily <- data.frame(dmvalue[,-c(1:2)])
    names(dmDaily) <- colnames(dm.data)
    rownames(dmDaily) <- unique(date)
  }

  return(dmDaily)
}
