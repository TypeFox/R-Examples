#' Create average waveform plots for each subject in a single, multiplot window
#'
#' \code{mosaic} generates multiple plots in a single window. Plots are average waveforms
#'   for each subject. Each plot shows all conditions present in the data frame.
#'
#' @param data A data frame in the format returned from \code{\link{load.data}}
#' @param electrodes A single value or concatenation of several values (to be averaged)
#'   indicating which electrodes to include in generating the plot. At this time, if the
#'   raw data files imported using \code{\link{load.data}}) do not have a header, you
#'   must include a capital "V" in front of the number and enclose each electrode in quotes.
#'   (For example, electrodes = "V78", or electrodes = c("V78", "V76").)
#' @param cols An integer defining the number of desired columns of plots. The default is 5.
#' @param rows An integer defining the number of desired rows of plots. The default is 5.
#'
#' @details The default values for columns and rows (i.e., 5 and 5) and higher are best suited to
#'   the graphical parameters specified in the code. At this time, graphical parameters (e.g.,
#'   tick marks and labels) do not scale with the number of rows and colums. As this feature is
#'   not intended to produce manuscript-ready plots, the option to explore your data with any
#'   number of rows and columns is available, but fewer rows and columns will yield crowded plots.
#'
#'   Single electrodes can be passed to the package functions, or several electrodes can
#'   be provided (i.e., when using dense arrays) and those electrodes will be averaged
#'   together as a single electrode.
#'
#' @return A single window containing multiple plots (1 per subject)
#'
#' @examples
#' # Inspect average ERP waveforms for each subject
#' mosaic(ERPdata, electrodes = "V78", cols = 5, rows = 5)
#'
#' @author Travis Moore

# function that plots the individual, average waveforms
# in separate windows
mosaic <- function(data, electrodes, cols = 5, rows = 5) {
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
  # calls the cluster.seg function
  cluster <- .cluster.seg(data, electrodes)
  # calls the avg.sub function
  avgsub <- .avg.subs(data, electrodes, window, cluster, Time.range, trial.types)
  # extracts grand mean data
  means.cond.sub <- .ind.by.cond(data, electrodes, window, Time.range,
                                avgsub, trial.types, Stimulus, num.subs, num.conditions)
  m = as.data.frame(means.cond.sub)

opar <- par(no.readonly = TRUE)
on.exit(par(opar))
  par(mfrow=c(rows, cols), mai=c(.1 , .1, .2, .1), mgp = c(-3,-1.5,0))
  for (i in 1:num.subs) {
    plot(unlist(m[i]) ~ Time.range,
         typ = "l",
         lwd = 1,
         col = 2,
         ylim = c(min(m),max(m)),
         main = sub.IDs[i],
         tck = 0.05,
         xlab = "",
         ylab = "")
    box()
    count <- i  # keeps track of how many times the main 'for' loop has run
    # (i.e., which subject number)
    counter <- 3
    for (j in c(seq(count, ncol(m), num.subs))) {
      counter <- counter + 1
      try(  # this catches a known error in creating the individual plots.
        # try() ~ on.error.resume.next from VB
        lines(unlist(m[ (j + num.subs)]) ~ Time.range,
              typ = "l",
              lwd = 1,
              col = counter),
        silent = TRUE)  # silent=TRUE with try() suppresses the error message(s)
    }  # close nested 'for' loop
  } # close 'for' loop
}  #close plot.ind()
