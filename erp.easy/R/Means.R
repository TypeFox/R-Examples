#' Calculate grand average and individual mean amplitudes and standard deviations
#'
#' \code{m.measures} calculates mean amplitude and standard deviation for each condition
#'   in the data frame, for the specified time window. Values are calculated based on grand
#'   average waveforms, as well as for each individual subject. Values are based on the electrode,
#'   or electrode cluster for dense arrays, provided in \code{electrodes}.
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
#'
#' @details Single electrodes can be passed to the package functions, or several electrodes can
#'   be provided (i.e., when using dense arrays) and those electrodes will be averaged
#'   together as a single electrode.
#'
#' @return A data frame with columns labeled:
#' \itemize{
#'     \item Subject
#'     \item Trial Type
#'     \item Standard Deviation
#'     \item Mean Amplitude
#' }
#'
#' @examples
#' # Calculate mean amplitude and standard deviation
#' m.measures(ERPdata, electrodes = "V78", window = c(1000, 1500))
#'
#' @author Travis Moore

m.measures <- function(data, electrodes, window) {
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
  win1 <- window[1]
  win2 <- window[2]
  # calls the cluster function
  cluster <- .cluster.seg(data, electrodes)
  # calls the avg.sub function
  avgsub <- .avg.subs(data, electrodes, window, cluster, Time.range, trial.types)
  # extracts grand mean data
  grand.avg <- .grand.average(data, electrodes, window, Time.range,
                              avgsub, trial.types, Stimulus, num.subs, num.conditions)

  # vectors to hold mean measures
# ------------- starts with grand average measures
means.ga <- vector("list", num.conditions)  # start with grand average measures
for (t in 1:num.conditions) {
  means.ga[[t]] <- lapply(trial.types[t], grand.avg, Stimulus, Time.range,
                          win1, win2, FUN = .get.ga.mamps)
}
sd.ga <- vector("list", num.conditions)
for (t in 1:num.conditions) {
  sd.ga[[t]] = lapply(trial.types[t], grand.avg, Stimulus, Time.range,
                      win1, win2, FUN = .get.ga.msds)
}
# get the peak values to use them in order to find the subject.  Cannot search raw value
# data frame for subject using means
peaks.ga <- vector("list", num.conditions)
for (t in 1:num.conditions) {
  peaks.ga[[t]] = lapply(trial.types[t], grand.avg, Stimulus, Time.range,
                         win1, win2, num.pts = 10, FUN = .get.ga.pamps)
}
unpacked.means.ga <- unlist(means.ga)
unpacked.sd.ga <- unlist(sd.ga)
unpacked.peaks.ga <- unlist(peaks.ga)
# comes afterward because it uses unpacked.peaks.ga
means.ga.cond = lapply(unpacked.peaks.ga, grand.avg, Stimulus, Time.range, win1, win2,
                       FUN = .get.peak.ga.cond)
unpacked.mean.ga.cond <- unlist(means.ga.cond)
grand.means <- setNames(data.frame("Grand Avg", unpacked.mean.ga.cond, unpacked.sd.ga,
                                   unpacked.means.ga), c("Subject", "Trial Type",
                                                      "Standard Dev", "Mean Amplitude"))

# ------------- begin individual measures
means.ind <- vector("list", num.conditions)
for (i in 1:num.conditions) {
  means.ind[[i]] = lapply(sub.IDs, trial.types[i], avgsub, Stimulus.ind, Time.range,
                          win1, win2, FUN = .get.mean.amps)
}
sd.mean <- vector("list", num.conditions)
for (t in 1:num.conditions) {
  sd.mean[[t]] <- lapply(sub.IDs, trial.types[t], avgsub, Stimulus.ind, Time.range,
                         win1, win2, FUN = .get.mean.msds)
}

peaks = vector("list", num.conditions)	# begin individual measures
for (i in 1:num.conditions) {
  peaks[[i]] = lapply(sub.IDs, trial.types[i], avgsub, Stimulus.ind, Time.range,
                      win1, win2, num.pts = 10, FUN = .get.peak.amps)
}
unpacked.peak.amp <- unlist(peaks)
unpacked.mean.amp <- unlist(means.ind)
unpacked.mean.sd <- unlist(sd.mean)

peak.sub = vector("list", length(unpacked.mean.amp))
peaks.sub = lapply(unpacked.peak.amp, avgsub, Stimulus.ind, Time.range,
                   win1, win2, FUN = .get.peak.sub)
unpacked.peak.sub <- unlist(peaks.sub)

peak.cond = vector("list", length(unpacked.peak.amp))
peaks.cond = lapply(unpacked.peak.amp, avgsub, Stimulus.ind, Time.range,
                    win1, win2, FUN = .get.peak.cond)
unpacked.peak.cond <- unlist(peaks.cond)

ind.means <- setNames(data.frame(unpacked.peak.sub, unpacked.peak.cond,
                                 unpacked.mean.sd, unpacked.mean.amp),
                      c("Subject", "Trial Type", "Standard Dev", "Mean Amplitude"))
mean.measures <- rbind(grand.means, ind.means)

grandaverage(data, electrodes, window)

return(mean.measures)
} # Close main function
