#' Define stem cycles and calculate statistics for all cyclic phases
#'
#' @description The function defines stem cycles from output of \code{\link{phase_def}} and calculates statistics for complete cycles as well as for the phases of contraction, expansion and stem-radius increment.
#'
#' @usage cycle_stats(dm.gpf, dm.phase, sensor = 1, smooth.param = 1)
#'
#' @param dm.gpf a \code{data.frame} with either gap-free or gap-filled dendrometer series as produced by \code{\link{fill_gaps}}.
#' @param dm.phase a \code{data.frame} with numbers indicating the different stem-cyclic phases. Output of \code{\link{phase_def}}.
#' @param sensor a \code{numeric} specifying the sensor to be used in the function. Defaults to 1 (first column in both \code{data.frames}).
#' @param smooth.param a \code{numeric} specifying the degree of smoothing. Defaults to 1 (no smoothing).
#'
#' @details The function uses the output of \code{\link{phase_def}} to define stem cycles and to calculate statistics for all cyclic phases. These statistics include the timing and duration of each phase, as well as information on stem-size changes. The function works for single dendrometer series, which are defined by the argument \code{sensor}.
#'
#' The function includes a smoothing option (argument \code{smooth.param}) particularly for noisy datasets in which outliers may under- or overestimate the minimum and maximum stem size within phases and stem cycles. By default, no smoothing is performed.
#'
#' @return The function returns a \code{list} with:
#' \itemize{\item{a \code{data.frame} named \code{cycleStats} containing the following summary statistics:}}
#' \item{dmID}{dendrometer ID.}
#' \item{cycle}{cycle number.}
#' \item{phase}{cyclic phase (1: contraction, 2: expansion, 3: stem-radius increment, 4: full cycle).}
#' \item{begin}{timestamp indicating the beginning of each phase.}
#' \item{end}{timestamp indicating the end of each phase.}
#' \item{duration_h}{phase duration in hours.}
#' \item{duration_m}{phase duration in minutes.}
#' \item{magnitude}{magnitude of stem-size changes in each phase.}
#' \item{min}{minimum stem size within each phase.}
#' \item{max}{maximum stem size within each phase.}
#'
#' \itemize{\item{a \code{data.frame} named \code{cycle.df} containing, for all individual records, the following columns:}}
#' \item{dmID}{dendrometer ID.}
#' \item{cycle}{cycle number.}
#' \item{phase}{cyclic phase (1: contraction, 2: expansion, 3: stem-radius increment, 4: full cycle).}
#'
#' @author Olivier Bouriaud, Ernst van der Maaten and Marieke van der Maaten-Theunissen.
#'
#' @examples
#' data(dmCD)
#' dm.phase <- phase_def(dmCD)
#' dm.stats <- cycle_stats(dmCD, dm.phase)
#'
#' @import stats
#'
#' @export cycle_stats
#'
cycle_stats <- function(dm.gpf, dm.phase, sensor = 1, smooth.param = 1)
{
  nm1 <- deparse(substitute(dm.gpf))
  nm2 <- deparse(substitute(dm.phase))

  if(!is.dendro(dm.gpf)) {
    stop(paste("'", nm1, "' is not in the required format", sep = ""))
  }
  if(!is.dendro(dm.phase)) {
    stop(paste("'", nm2, "' is not in the required format", sep = ""))
  }
  if(nrow(dm.gpf) != nrow(dm.phase)) {
    stop(paste("'", nm1, "' and '", nm2,"' have different numbers of rows", sep = ""))
  }
  if(ncol(dm.gpf) != ncol(dm.phase)) {
    stop(paste("'", nm1, "' and '", nm2,"' have different numbers of columns", sep = ""))
  }
  if(ncol(dm.gpf) == 1 && sensor > ncol(dm.gpf)) {
    stop("'sensor' should be 1")
  }
  if(ncol(dm.gpf) > 1 && sensor > ncol(dm.gpf)) {
    stop(paste("'sensor' should be between", 1, "and", ncol(dm.gpf), sep = " "))
  }
  if(dendro.resolution(dm.phase)!=dendro.resolution(dm.gpf)) {
    stop(paste("'", nm1, "' and '", nm2, "' do not have the same temporal resolution", sep = ""))
  }

  # Optional smoothing
  if(smooth.param > 1) {
    dm.gpf[, sensor] <- ave(dm.gpf[, sensor],
                            FUN = function(x) rollmean(x, smooth.param, align = "right", fill = NA, na.rm = TRUE))
  }

  dm.gpj <- cbind(dm.phase[, sensor, drop = FALSE], dm.gpf[, sensor, drop = FALSE])
  dm.gpj <- na.omit(dm.gpj)
  date <- strptime(row.names(dm.gpj), "%Y-%m-%d %H:%M:%S")

  # Cycle definition
  cycle <- vector(length = dim(dm.gpj)[1])
  k <- 1
  for(i in 2:(dim(dm.gpj)[1])) {
    if(is.na(dm.gpj[i,1]) == FALSE) {
      if(dm.gpj[i,1] == 1 & dm.gpj[i-1,1]>1) {
        k <- k+1
      }
    }
    cycle[i] <- k
  }

  # Calculate cycle stats
  frq.sec <- dendro.resolution(dm.gpf, unts = "secs")

  phase <- dm.gpj[,1]
  date.begin <- aggregate(date, list(phase, cycle), FUN = min)$x - frq.sec
  date.end <- aggregate(date, list(phase, cycle), FUN = max)$x
  duration.h <- difftime(date.end, date.begin, units = "hours")
  duration.m <- difftime(date.end, date.begin, units = "mins")
  cycle.min <- aggregate(dm.gpj[,2], list(phase, cycle), FUN = min)$x
  cycle.max <- aggregate(dm.gpj[,2], list(phase, cycle), FUN = max)$x
  magnitude <- cycle.max - cycle.min
  cycle.nbr <- unique(data.frame(cycle = cycle, phase = phase))

  dm.stats <- data.frame(dmID = rep(names(dm.gpj)[2], length(date.begin)),
                         cycle = cycle.nbr$cycle,
                         phase = cycle.nbr$phase,
                         begin = date.begin,
                         end = date.end,
                         duration_h = duration.h,
                         duration_m = duration.m,
                         magnitude = magnitude,
                         min = cycle.min,
                         max = cycle.max)

  cycle4 <- unique(dm.stats$cycle)
  phase4 <- rep(4, l = length(cycle4))
  date.begin4 <- aggregate(date, list(cycle, cycle), FUN = min)$x - frq.sec
  date.end4 <- aggregate(date, list(cycle, cycle), FUN = max)$x
  duration.h4 <- difftime(date.end4, date.begin4, units = "hours")
  duration.m4 <- difftime(date.end4, date.begin4, units = "mins")
  cycle.min4 <- aggregate(dm.gpj[,2], list(cycle, cycle), FUN = min)$x
  cycle.max4 <- aggregate(dm.gpj[,2], list(cycle, cycle), FUN = max)$x
  magnitude4 <- cycle.max4 - cycle.min4

  dm.stats4 <- data.frame(dmID = rep(names(dm.gpj)[2], length(date.begin4)),
                          cycle = cycle4,
                          phase = phase4,
                          begin = date.begin4,
                          end = date.end4,
                          duration_h = duration.h4,
                          duration_m = duration.m4,
                          magnitude = magnitude4,
                          min = cycle.min4,
                          max = cycle.max4)

  dm.stats <- rbind(dm.stats, dm.stats4)
  dm.stats <- dm.stats[dm.stats$cycle>0, ]
  dm.stats <- dm.stats[order(dm.stats$cycle, dm.stats$phase),]
  rownames(dm.stats) <- seq(1, nrow(dm.stats), 1)

  cycle.df <- data.frame(dmID = rep(names(dm.gpj)[2], length(cycle)),
                         cycle = cycle, phase = phase)
  rownames(cycle.df) <- date
  output <- list(cycleStats = dm.stats, cycle.df = cycle.df)
  return(output)
}
