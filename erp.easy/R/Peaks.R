#' Calculate grand average and individual peak amplitude and latency
#'
#' \code{p.measures} calculates simple peak amplitude and latency for each condition
#'   in the data frame. Values are calculated for grand average waveforms, as well as
#'   for each individual subject. Values are based on the electrode, or electrode cluster for
#'   dense arrays, provided in \code{electrodes}. This function will identify the largest
#'   deviation from 0, whether positive or negative.
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
#' @param num.pts The number of bins to check for local peak measures. If no local peaks are
#'   found, the simple peak will be returned. To force the simple peak, set \code{num.pts}
#'   to 0.
#'
#' @details At this time there is no way to specify a negative or positive peak.
#'   \code{p.measures} simply returns the largest absolute value deviation from zero.
#'
#' @return A data frame with columns labeled:
#' \itemize{
#'     \item Subject
#'     \item Trial Type
#'     \item Peak Latency
#'     \item Peak Amplitude
#' }
#'
#' @examples
#' # Calculate peak latency and amplitude
#' p.measures(ERPdata, electrodes = "V78", window = c(1000, 1500))
#'
#' @author Travis Moore


p.measures <- function(data, electrodes, window, num.pts = 10) {
  data.fun <- data
  num.subs <- length(levels(data$Subject))
  sub.IDs <- levels(data$Subject)
  num.conditions <- length(levels(data$Stimulus))
  trial.types <- levels(data$Stimulus)
  time.points <- (length(data$Time)/num.subs)/num.conditions  # an integer
  # of the number of time points for one stimulus type for one subject
  Time.range <- data$Time[1:time.points]

# only use this for grand average measures - for alternating stimuli blocks (old way)
  stim.list <- vector("list")
  for (i in 1:length(levels(data$Stimulus))) {
    stim.list[[i]] <- c(rep(levels(data$Stimulus)[i], time.points))
  }
  Stimulus <- unlist(stim.list)

# only use this for independant measures
  Stimulus.ind <- data$Stimulus

  win1 <- window[1]
  win2 <- window[2]
  # calls the cluster function
  cluster <- .cluster.seg(data, electrodes)
  # calls the avg.sub function
  avgsub <- .avg.subs(data, electrodes, window, cluster, Time.range, trial.types)
  # extracts grand mean data
  grand.avg <- .grand.average(data, electrodes, window, Time.range,
                                  avgsub, trial.types, Stimulus, num.subs, num.conditions)

    # vectors to hold peak measures
# ------------- starts with grand average measures
  peaks.ga <- vector("list", num.conditions)
    # these variables would not be necessary if peaks.ga collected the whole row
    # but it was easier to program each piece of info for the data frame
    for (t in 1:num.conditions) {
      peaks.ga[[t]] = lapply(trial.types[t], grand.avg, Stimulus, Time.range,
                             win1, win2, num.pts, FUN = .get.ga.pamps)
    }
  unpacked.peaks.ga <- unlist(peaks.ga)
  peak.ga.lat = vector("list", length(unpacked.peaks.ga))
  peak.ga.cond = vector("list", length(unpacked.peaks.ga))
  peaks.ga.cond = lapply(unpacked.peaks.ga, grand.avg, Stimulus, Time.range, win1, win2,
                         FUN = .get.peak.ga.cond)
  peaks.ga.lat = lapply(unpacked.peaks.ga, grand.avg, Stimulus, Time.range, win1, win2,
                        FUN = .get.peak.ga.latency)
  unpacked.peak.ga.cond <- unlist(peaks.ga.cond)
  unpacked.peak.ga.lat <- unlist(peaks.ga.lat)
  grand.peaks = setNames(data.frame("Grand Avg", unpacked.peak.ga.cond,
                                    unpacked.peak.ga.lat, unpacked.peaks.ga),
                         c("Subject", "Trial Type", "Peak Latency", "Peak Amplitude"))

# ------------- begin individual measures

  # this is likely the culprit.  Make sure the lapply function is actually cycling through every subject
    peaks = vector("list",num.conditions)	# begin individual measures
  for (i in 1:num.conditions) {
    peaks[[i]] = lapply(sub.IDs, trial.types[i], avgsub, Stimulus.ind, Time.range,
                        win1, win2, num.pts, FUN = .get.peak.amps)


    #peaks[[i]] = lapply(sub.IDs, trial.types[i], avgsub, Stimulus.ind, Time.range,
    #                    win1, win2, num.pts, FUN = .get.peak.amps)
  } # close for loop

  unpacked.peak.amp <- unlist(peaks)
  peak.sub = vector("list", length(unpacked.peak.amp))
  peak.cond = vector("list", length(unpacked.peak.amp))
  peak.lat = vector("list", length(unpacked.peak.amp))
  peaks.sub = lapply(unpacked.peak.amp, avgsub, Stimulus.ind, Time.range,
                     win1, win2, FUN = .get.peak.sub)
  peaks.cond = lapply(unpacked.peak.amp, avgsub, Stimulus.ind, Time.range,
                      win1, win2, FUN = .get.peak.cond)
  peaks.lat = lapply(unpacked.peak.amp, avgsub, Stimulus.ind, Time.range,
                     win1, win2, FUN = .get.peak.latency)
    unpacked.peak.lat <- unlist(peaks.lat)
    unpacked.peak.sub <- unlist(peaks.sub)
    unpacked.peak.cond <- unlist(peaks.cond)

  peaks.measures = setNames(data.frame(unpacked.peak.sub, unpacked.peak.cond,
                                       unpacked.peak.lat, unpacked.peak.amp),
                            c("Subject", "Trial Type", "Peak Latency",
                              "Peak Amplitude"))	# individual measures df

  grandaverage(data, electrodes, window)
  tryCatch(abline(v = unpacked.peak.ga.lat, lty = 6, col = c(2, 4, seq(5, num.conditions+3, 1))),
           error = function(e) {abline(v = unpacked.peak.ga.lat, lty = 6, col = 2)}, silent = TRUE)

# binds both peak measures data frames and makes it available outside the function
  peak.measures <- rbind(grand.peaks, peaks.measures)
  return(peak.measures)

}  # Function closing bracket
