#' Segmenting climate and environmental data
#'
#' @description The function calculates means or sums, or extracts minimum or maximum values of environmental parameters for stem-cyclic phases as defined using \code{\link{cycle_stats}}.
#'
#' @usage climate_seg(env.data, dm.stats, value = c("mean", "min",
#'             "max", "sum"))
#'
#' @param env.data a \code{data.frame} with with a timestamp (\code{\%Y-\%m-\%d \%H:\%M:\%S} format) as row names, and a certain climate parameter (e.g., temperature or precipitation) in columns.
#' @param dm.stats a \code{list} as produced by \code{\link{cycle_stats}}.
#' @param value a \code{character} string of \code{"mean"}, \code{"min"}, \code{"max"} or \code{"sum"}, specifying whether means (e.g., for temperature) or sums (e.g., for precipitation) should be calculated, or minimum or maximum values should be extracted. Defaults to \code{"mean"}. Argument matching is performed.
#'
#' @details The function segments environmental parameters according to the stem-cyclic phases as defined using \code{\link{cycle_stats}}. Means, sums, and minimum and maximum values can be calculated or extracted.
#'
#' \code{env.data} should cover at least the same period as the dendrometer data used to define the cyclic phases, and should have the same (or a higher) temporal resolution.
#'
#' @return The function returns a \code{data.frame} with segmented environmental data. The \code{data.frame} contains the following columns:
#' \item{dmID}{dendrometer ID.}
#' \item{cycle}{cycle number.}
#' \item{phase}{cyclic phase (1: contraction, 2: expansion, 3: stem-radius increment, 4: full cycle).}
#' \item{begin}{timestamp indicating the beginning of each phase.}
#' \item{end}{timestamp indicating the end of each phase.}
#' \item{...}{columns with segmented environmental data (mean, min, max or sum).}
#'
#' @examples
#' \dontrun{
#'
#' data(dmED)
#' dm.gpf <- fill_gaps(dmED)
#' dm.phase <- phase_def(dm.gpf)
#' dm.stats <- cycle_stats(dm.gpf, dm.phase)
#' data(envED)
#' clim.phase <- climate_seg(envED, dm.stats, value = "mean")
#' }
#'
#' @import stats
#'
#' @export climate_seg
#'
climate_seg <- function(env.data, dm.stats, value = c("mean", "min", "max", "sum"))
{
  cycle.df <- dm.stats[[2]]
  nm1  <- deparse(substitute(env.data))
  nm2 <- deparse(substitute(dm.stats))
  if(!is.dendro(env.data)) {
    stop(paste("'", nm1, "' is not in the required format", sep = ""))
  }

  date.env.data <- strptime(row.names(env.data), "%Y-%m-%d %H:%M:%S")
  date.cycle <- strptime(rownames(cycle.df), "%Y-%m-%d %H:%M:%S")
  if(!all(date.cycle %in% date.env.data)) {
    stop(paste("the environmental data does not contain all timestamps of '", nm2, "'", sep = ""))
  }

  env.data2 <- env.data
  env.data2$date <- date.env.data
  cycle.df$date <- date.cycle
  comb.env.dm <- merge(env.data2, cycle.df, by = "date", all.x = TRUE)

  value <- match.arg(value, c("mean", "min", "max",  "sum"))
  ph123 <- aggregate(env.data, list(phase = comb.env.dm$phase, cycle = comb.env.dm$cycle), FUN = value)
  ph4 <- aggregate(env.data, list(cycle = comb.env.dm$cycle, cycle = comb.env.dm$cycle), FUN = value)

  cycleStats <- dm.stats[[1]]
  dm.stats123 <- cycleStats[which(cycleStats$phase == 1 | cycleStats$phase == 2 | cycleStats$phase == 3),]
  climPhase123 <- data.frame(dm.stats123[,1:5], ph123[ph123$cycle>0,-c(1:2)])
  climPhase4 <- data.frame(cycleStats[which(cycleStats$phase == 4),1:5], ph4[ph4$cycle>0,-c(1:2)])

  climPhase <- rbind(climPhase123, climPhase4)
  climPhase <- climPhase[order(climPhase$cycle, climPhase$phase),]
}
