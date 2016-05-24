#' Generate a butterfly plot for a specified condition
#'
#' \code{butterfly} plots all individual waveforms for the condition specified by the
#'   \code{stim} argument (i.e., a butterfly plot). The grand average waveform is also plotted,
#'   using a red line.
#'
#' @param data A data frame in the format returned from \code{\link{load.data}}
#' @param electrodes A single value or concatenation of several values (to be averaged)
#'   indicating which electrodes to include in generating the plot. At this time, if the
#'   raw data files imported using \code{\link{load.data}}) do not have a header, you
#'   must include a capital "V" in front of the number and enclose each electrode in quotes.
#'   (For example, electrodes = "V78", or electrodes = c("V78", "V76").)
#' @param stim An integer specifying which condition to plot. Conditions are numbered in the
#'   order in which they are imported with \code{load.data}.  The plot title will specify
#'   the name of the condition to avoid confusion.
#'
#' @details Single electrodes can be passed to the package functions, or several electrodes can
#'   be provided (i.e., when using dense arrays) and those electrodes will be averaged
#'   together as a single electrode.
#'
#' @return A single butterfly plot for the condition specified with \code{stim}.
#'
#' @examples
#' butterfly(ERPdata, electrodes = "V78", stim = 1)
#'
#' @author Travis Moore


# function that plots all the individual waveforms
# for the specified condition, overlaid also plots the
# grand average for that condition in red
butterfly <- function(data, electrodes, stim = 1) {
  data.fun <- data
  num.subs <- length(levels(data$Subject))
  sub.IDs <- levels(data$Subject)
  num.conditions <- length(levels(data$Stimulus))
  trial.types <- levels(data$Stimulus)
  time.points <- (length(data$Time)/num.subs)/num.conditions  # an integer
  # of the number of time points for one stimulus type for one subject
  Time.range <- data$Time[1:time.points]
  # only use stim.list for grand average measures - for alternating stimuli blocks (old way)
  stim.list <- vector("list")
  for (i in 1:length(levels(data$Stimulus))) {
    stim.list[[i]] <- c(rep(levels(data$Stimulus)[i], time.points))
  }
  Stimulus <- unlist(stim.list)
  # only use this for independant measures
  Stimulus.ind <- data$Stimulus
  # calls the cluster function
  cluster <- .cluster.seg(data, electrodes)
  # calls the avg.sub function
  avgsub <- .avg.subs(data, electrodes, window, cluster, Time.range, trial.types)
  # extracts grand mean data
  grand.avg <- .grand.average(data, electrodes, window, Time.range,
                              avgsub, trial.types, Stimulus, num.subs, num.conditions)
  means.cond.sub <- .ind.by.cond(data, electrodes, window, Time.range,
                                 avgsub, trial.types, Stimulus, num.subs, num.conditions)
  plot(unlist(means.cond.sub[[stim]][ , 1]) ~ Time.range,
       typ = "l", ylim = c(min(unlist(means.cond.sub[[stim]])),
                           max(unlist(means.cond.sub[[stim]]))),
       xlab = "Time (ms)",
       ylab = "Amplitude (microV)",
       main = paste("Individual waveforms for", trial.types[stim], seq = ""))

  if (num.subs > 1) {
    #for (i in 2:num.subs) {
    for (i in 2:num.subs) {
      lines(unlist(means.cond.sub[[stim]][ , i]) ~ Time.range, typ = "l")
    }

    with(grand.avg, {
  grand.line <- subset(grand.avg, Stimulus == trial.types[stim], select = Means)
  lines(unlist(grand.line) ~ Time.range,
        typ = "l",
        lwd = 2,
        col = "red")
    }
    )


  }
}  # close butterfly.plot()
