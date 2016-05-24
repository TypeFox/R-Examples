#' Plot the grand average waveform for all loaded conditions
#'
#' \code{grandaverage} plots the grand average waveform for each condition present in the
#'   data frame you provide.  A color-coded and labeled legend is generated with the plot
#'   for ease of identification of each condition.
#'
#' @param data A data frame in the format returned from \code{\link{load.data}}
#' @param electrodes A single value or concatenation of several values (to be averaged)
#'   indicating which electrodes to include in generating the plot. At this time, if the
#'   raw data files imported using \code{\link{load.data}}) do not have a header, you
#'   must include a capital "V" in front of the number and enclose each electrode in quotes.
#'   (For example, electrodes = "V78", or electrodes = c("V78", "V76").)
#' @param window The beginning and end points of a time window of interest; this is different
#'   from the beginning and ending times \code{epoch.st} and \code{epoch.end} defined in
#'   \code{\link{load.data}} (you only need to define the epoch once upon importing the data).
#'   The purpose of the \code{window} argument in this function is merely to highlight a
#'   window of interest; its default value is NULL.
#'
#' @details \code{grandaverage} will return a plot of the grand average waveform for each
#'   condition present in the data frame you provide.  For ease of use, colors are
#'   automatically assigned. The legend displays the value provided in the \code{condition}
#'   argument of \code{\link{load.data}}.
#'
#'   Single electrodes can be passed to the package functions, or several electrodes
#'   can be provided (i.e., when using dense arrays) and those electrodes will be
#'   averaged together as a single electrode.
#'
#' @return A single plot of grand average waveforms for each condition. Includes a
#'   color-coded and labeled legend.
#'
#' @examples
#' # Create a plot of the grand average waveforms for each imported condition
#' grandaverage(ERPdata, electrodes = "V78", window = c(1000, 1500))
#'
#' @author Travis Moore

  # plots the grand average of all loaded conditions
grandaverage <- function(data, electrodes, window = NULL) {
  data.fun <- data
  num.subs <- length(levels(data$Subject))
  sub.IDs <- levels(data$Subject)
  num.conditions <- length(levels(data$Stimulus))
  trial.types <- levels(data$Stimulus)
  stim.block <- length(data$Time)/num.subs
  Stimulus <- data$Stimulus[1:stim.block]
  time.points <- (length(data$Time)/num.subs)/num.conditions  # an integer
    # of the number of time points for one stimulus type for one subject
  Time.range <- data$Time[1:time.points]
    if (!is.null(window)) {
  win1 <- window[1]
  win2 <- window[2]
    } else {
  # nothing
    }
    # calls the cluster function
  cluster <- .cluster.seg(data, electrodes)
    # calls the avg.sub function
  avgsub <- .avg.subs(data, electrodes, window, cluster, Time.range, trial.types)
    # extracts grand mean data
  means.cond.sub <- .ind.by.cond(data, electrodes, window, Time.range,
                            avgsub, trial.types, Stimulus, num.subs, num.conditions)
    # does the actual plotting
  lower <- min(rowMeans(plyr::ldply(means.cond.sub)))
  upper <- max(rowMeans(plyr::ldply(means.cond.sub)))
  plot(rowMeans(as.data.frame(means.cond.sub[1])) ~ Time.range, typ = "l", lwd = 3,
       main = "Grand Average", xlab = "Time in milliseconds", ylab = "Amplitude in microvolts",
       col = 2, ylim = c(lower, upper))



  if (!is.null(window)) {
    rect(win1, lower, win2, upper, lwd = 2, col=rgb(0.8, 0.8, 0.8, 0.3))

    #legend("topright", inset = 0.05, title = "Trial Types", lwd = 3, trial.types,
    #      col = 2, bg = "white")

  } else {
    # nothing
  }




  legend("topright", inset = 0.05, title = "Trial Types", lwd = 3, trial.types,
         col = 2, bg = "white")

    if (num.conditions >= 2) {
      counter <- 3 # counter that increments the color of the lines()
      for (k in 2:num.conditions) {
        counter <- counter + 1
        lines(rowMeans(as.data.frame(means.cond.sub[k])) ~ Time.range, typ = "l",
              lwd = 3, col = counter)

        try(legend("topright", inset = 0.05, title = "Trial Types", lwd = 3, trial.types,
                   col = c(2, 4, seq(5, num.conditions+3, 1)), bg = "white")
            , silent = TRUE)  # suppresses a known error when there is only 1 condition
      }
    } else {
      # nothing
    }

  # activate this line to return a data frame of the raw amplitudes of the grand average
  return(rowMeans(as.data.frame(means.cond.sub)))


}  # CLOSE MAIN FUNCTION
