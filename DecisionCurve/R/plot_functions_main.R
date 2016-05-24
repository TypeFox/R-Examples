#'Plot the net benefit curves from a decision_curve object or many decision_curve objects
#'
#' @param x 'decision_curve' object to plot or a list of 'decision_curve' objects. Assumes output from function 'decision_curve'
#' @param curve.names vector of names to use when plotting legends.
#' @param cost.benefit.axis logical (default TRUE) indicating whether to print an additional x-axis showing relative cost:benefit ratios in addition to risk thresholds.
#' @param n.cost.benefits number of cost:benefit ratios to print if cost.benefit.axis = TRUE (default n.cost.benefit = 6).
#' @param cost.benefits Character vector of the form c("c1:b1", "c2:b2", ..., "cn:bn") with integers ci, bi corresponding to specific cost:benefit ratios to print. Default allows the function to calculate these automatically.
#' @param standardize logical (default TRUE) indicating whether to use the standardized net benefit (NB/disease prevalence) or not.
#' @param confidence.intervals logical indicating whether to plot confidence intervals.
#' @param col vector of color names to be used in plotting corresponding to the 'predictors' given. Default colors will be chosen from rainbow(..., v = .8). See details for more information on plot parameters.
#' @param lty vector of linetypes.
#' @param lwd vector of linewidths.
#' @param xlim vector giving c(min, max) of x-axis. Defaults to c(min(thresholds), max(thresholds)).
#' @param ylim vector giving c(min, max) of y-axis.
#' @param xlab label of main x-axis.
#' @param ylab label of y-axis.
#' @param cost.benefit.xlab label of cost:benefit ratio axis.
#' @param legend.position character vector giving position of legend. Options are "topright" (default), "right", "bottomright", "bottom", "bottomleft", "left", "topleft", "top", or "none".
#' @param ... other options directly send to plot()
#' @details When k decision_curve objects are input, the first k elements of col, lty, lwd ... correspond to the curves provided. The next two elements (..., k+1, k+2) correspond to the attributes of the 'all' and 'none' curves. See below for an example.
#' @examples
#'data(dcaData)
#'set.seed(123)
#'baseline.model <- decision_curve(Cancer~Age + Female + Smokes,
#'                                 data = dcaData,
#'                                 thresholds = seq(0, .4, by = .005),
#'                                 bootstraps = 10)
#'
#'#plot using the defaults
#'plot_decision_curve(baseline.model,  curve.names = "baseline model")
#'
#'set.seed(123)
#'full.model <- decision_curve(Cancer~Age + Female + Smokes + Marker1 + Marker2,
#'                             data = dcaData,
#'                             thresholds = seq(0, .4, by = .005),
#'                             bootstraps = 10)
#'
#'# for lwd, the first two positions correspond to the decision curves, then 'all' and 'none'
#'plot_decision_curve( list(baseline.model, full.model),
#'                     curve.names = c("Baseline model", "Full model"),
#'                     col = c("blue", "red"),
#'                     lty = c(1,2),
#'                     lwd = c(3,2, 2, 1),
#'                     legend.position = "bottomright")
#'
#No confidence intervals, cost:benefit ratio axis, or legend
#'
#'plot_decision_curve( list(baseline.model, full.model),
#'                     curve.names = c("Baseline model", "Full model"),
#'                     col = c("blue", "red"),
#'                     confidence.intervals = FALSE,  #remove confidence intervals
#'                     cost.benefit.axis = FALSE, #remove cost benefit axis
#'                     legend.position = "none") #remove the legend
#'
#'#Set specific cost:benefit ratios.
#'
#'plot_decision_curve( list(baseline.model, full.model),
#'                     curve.names = c("Baseline model", "Full model"),
#'                     col = c("blue", "red"),
#'                     cost.benefits = c("1:1000", "1:4", "1:9", "2:3", "1:3"),
#'                     legend.position = "bottomright")
#'
#'#Plot net benefit instead of standardize net benefit.
#'
#'plot_decision_curve( list(baseline.model, full.model),
#'                     curve.names = c("Baseline model", "Full model"),
#'                     col = c("blue", "red"),
#'                     ylim = c(-0.05, 0.15), #set ylim
#'                     lty = c(2,1),
#'                     standardize = FALSE, #plot Net benefit instead of standardized net benefit
#'                    legend.position = "topright")
#'
#'
#' @export

plot_decision_curve <- function(x, curve.names,
                               cost.benefit.axis = TRUE,
                               n.cost.benefits = 6,
                               cost.benefits,
                               standardize = TRUE,
                               confidence.intervals,
                               col,
                               lty, lwd = 2,
                               xlim, ylim,
                               xlab = "Risk Threshold", ylab,
                               cost.benefit.xlab = "Cost:Benefit Ratio",
                               legend.position = c("topright", "right", "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "none"),
                               ...){

  legend.position <- match.arg(legend.position)

  if(missing(curve.names)) curve.names  <- NA
  if(missing(confidence.intervals)) confidence.intervals <- NA

  prepData <- preparePlotData(x = x,
                             curve.names = curve.names,
                             confidence.intervals = confidence.intervals)

  predictors <- prepData$predictors
  dc.data <- prepData$dc.data
  confidence.intervals <- prepData$confidence.intervals
  rm(prepData)

  #set some defaults if needed
  if(missing(xlim)) xlim = range(dc.data$thresholds)

  if(missing(lty)) lty = rep(1, length(predictors) + 2)
  if(length(lty) ==1) lty = rep(lty, length(predictors) + 2)
  if(length(lty) == length(predictors)) lty = c(lty, 1, 1)

  if(missing(col)) col  = c(rainbow(length(predictors), v = .8), "grey66", "black")
  if(length(col) == length(predictors)) col <- c(col, "grey66", "black")

  if(missing(lwd)) lwd = 2
  if(length(lwd) ==1) lwd <- rep(lwd, length(predictors))
  if(length(lwd) == length(predictors)) lwd = c(lwd, 1, 1)

  if(missing(ylab)) ylab <- ifelse(standardize, "Standardized Net Benefit", "Net Benefit")

    if(missing(ylim)){

    if(standardize) ylim = c(-1, 1)
    else ylim = c(-0.05, 1.1*max(dc.data[["NB"]][is.finite(dc.data[["NB"]])]))

  }

  plot_generic(xx = dc.data,
               predictors = predictors,
               value = ifelse(standardize, "sNB", "NB"),
               plotNew = TRUE,
               standardize = standardize,
               confidence.intervals,
               cost.benefit.axis = cost.benefit.axis,
               cost.benefits = cost.benefits,
               n.cost.benefits = n.cost.benefits,
               cost.benefit.xlab = cost.benefit.xlab,
               xlab = xlab, ylab = ylab,
               col = col,
               lty = lty, lwd = lwd,
               xlim = xlim, ylim = ylim,
               legend.position = legend.position, ...)

}

#'Plot the components of a ROC curve by the high risk thresholds.
#' @description Plot the components of the ROC curve --the true positive rates and false positive rates-- by high risk thresholds.
#' @param x decision_curve object to plot. Assumes output from function 'decision_curve'
#' @param cost.benefit.axis logical (default TRUE) indicating whether to print an additional x-axis showing relative cost:benefit ratios in addition to risk thresholds.
#' @param n.cost.benefits number of cost:benefit ratios to print if cost.benefit.axis = TRUE (default n.cost.benefit = 6).
#' @param cost.benefits Character vector of the form c("c1:b1", "c2:b2", ..., "cn:bn") with integers ci, bi corresponding to specific cost:benefit ratios to print. Default allows the function to calculate these automatically.
#' @param confidence.intervals logical indicating whether to plot confidence intervals.
#' @param col vector of length two indicating the color for the true positive rates and false positive rates, respectively.
#' @param lty.fpr linetype for the false positive rate curve.
#' @param lty.tpr linetype for the true positive rate curve.
#' @param lwd vector of linewidths. The first element corresponds to the tpr and the second to the fpr.
#' @param xlim vector giving c(min, max) of x-axis. Defaults to c(min(thresholds), max(thresholds)).
#' @param ylim vector giving c(min, max) of y-axis.
#' @param xlab label of main x-axis.
#' @param ylab label of y-axis.
#' @param cost.benefit.xlab label of cost:benefit ratio axis.
#' @param legend.position character vector giving position of legend. Options are "topright" (default), "right", "bottomright", "bottom", "bottomleft", "left", "topleft", "top", or "none".
#' @param ... other options directly send to plot()
#'
#' @examples
#'data(dcaData)
#'set.seed(123)
#'baseline.model <- decision_curve(Cancer~Age + Female + Smokes,
#'                                 data = dcaData,
#'                                 thresholds = seq(0, .4, by = .001),
#'                                 bootstraps = 25) #should use more bootstrap replicates in practice!
#'
#'#plot using the defaults
#'plot_roc_components(baseline.model,  xlim = c(0, 0.4), col = c("black", "red"))
#'
#'
#' @export


plot_roc_components <- function(x,
                              cost.benefit.axis = TRUE,
                              n.cost.benefits = 6,
                              cost.benefits,
                              confidence.intervals,
                              col = "black",
                              lty.fpr = 2,
                              lty.tpr = 1,
                              lwd = 2,
                              xlim, ylim,
                              xlab = "Risk Threshold", ylab,
                              cost.benefit.xlab = "Cost:Benefit Ratio",
                              legend.position = c("topright", "right", "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "none"),
                              ...){

  if(class(x) != "decision_curve") stop("x must be an object of class 'decision_curve'-- plot_roc_components is only functional for one decision curve at a time.")

  legend.position <- match.arg(legend.position)

 # if(missing(curve.names)) curve.names  <- NA
  if(missing(confidence.intervals)) confidence.intervals <- NA

  prepData <- preparePlotData(x = x,
                              curve.names = NA,
                              confidence.intervals = confidence.intervals)

  predictors <- prepData$predictors
  dc.data <- prepData$dc.data
  confidence.intervals <- prepData$confidence.intervals
  rm(prepData)

  #set some defaults if needed
  if(missing(xlim)) xlim = range(dc.data$thresholds)

  if(length(col) ==1 ) col <- rep(col, 2)

  if(missing(lwd)) lwd = 2
  if(length(lwd) ==1) lwd <- rep(lwd, length(predictors))
  if(length(lwd) == length(predictors)) lwd = c(lwd, 1, 1)

  if(missing(ylab)) ylab <- "Probability"

  if(missing(ylim)) ylim = c(0, 1)


  plot_generic(xx = dc.data,
               predictors = predictors,
               value = "TPR",
               plotNew = TRUE,
               standardize = FALSE,
               confidence.intervals,
               cost.benefit.axis = cost.benefit.axis,
               cost.benefits = cost.benefits,
               n.cost.benefits = n.cost.benefits,
               cost.benefit.xlab = cost.benefit.xlab,
               xlab = xlab, ylab = ylab,
               col = col,
               lty = lty.tpr, lwd = lwd,
               xlim = xlim, ylim = ylim,
               legend.position = "none",
               ...)

  plot_generic(xx = dc.data,
               predictors = predictors,
               value = "FPR",
               plotNew = FALSE,
               standardize = FALSE,
               confidence.intervals,
               cost.benefit.axis = cost.benefit.axis,
               cost.benefits = cost.benefits,
               n.cost.benefits = n.cost.benefits,
               cost.benefit.xlab = cost.benefit.xlab,
               xlab = xlab, ylab = ylab,
               col = col,
               lty = lty.fpr, lwd = lwd,
               xlim = xlim, ylim = ylim,
               legend.position = legend.position,
               lty.fpr = lty.fpr,
               lty.tpr = lty.tpr,
               tpr.fpr.legend = TRUE,
               ...)



}

#'Plot the clinical impact curve from a DecisionCurve object.
#'
#' @description For a given population size, plot the number of subjects classified as high risk, and the number of subjects classified high risk with the outcome of interest at each high risk threshold.
#' @param x decision_curve object to plot. Assumes output from function 'decision_curve'
#' @param population.size Hypothetical population size (default 1000).
#' @param cost.benefit.axis logical (default TRUE) indicating whether to print an additional x-axis showing relative cost:benefit ratios in addition to risk thresholds.
#' @param n.cost.benefits number of cost:benefit ratios to print if cost.benefit.axis = TRUE (default n.cost.benefit = 6).
#' @param cost.benefits Character vector of the form c("c1:b1", "c2:b2", ..., "cn:bn") with integers ci, bi corresponding to specific cost:benefit ratios to print. Default allows the function to calculate these automatically.
#' @param confidence.intervals logical indicating whether to plot confidence intervals.
#' @param col vector of length two indicating the color for the number high risk and the second to the number high risk with outcome, respectively.
#' @param lty vector of linetypes. The first element corresponds to the number high risk and the second to the number high risk with outcome.
#' @param lwd vector of linewidths. The first element corresponds to the number high risk and the second to the number high risk with outcome.
#' @param xlim vector giving c(min, max) of x-axis. Defaults to c(min(thresholds), max(thresholds)).
#' @param ylim vector giving c(min, max) of y-axis.
#' @param xlab label of main x-axis.
#' @param ylab label of y-axis.
#' @param cost.benefit.xlab label of cost:benefit ratio axis.
#' @param legend.position character vector giving position of legend. Options are "topright" (default), "right", "bottomright", "bottom", "bottomleft", "left", "topleft", "top", or "none".
#' @param ... other options directly send to plot()
#'
#' @examples
#'#'data(dcaData)
#'set.seed(123)
#'baseline.model <- decision_curve(Cancer~Age + Female + Smokes,
#'                                 data = dcaData,
#'                                 thresholds = seq(0, .4, by = .001),
#'                                 bootstraps = 25) #should use more bootstrap replicates in practice!
#'
#'#plot the clinical impact
#'plot_clinical_impact(baseline.model, xlim = c(0, .4),
#'                     col = c("black", "blue"))
#'
#' @export

plot_clinical_impact <- function(x,
                                 population.size = 1000,
                            cost.benefit.axis = TRUE,
                            n.cost.benefits = 6,
                            cost.benefits,
                            confidence.intervals,
                            col = "black",
                            lty = 1,
                            lwd = 2,
                            xlim, ylim,
                            xlab = "Risk Threshold", ylab,
                            cost.benefit.xlab = "Cost:Benefit Ratio",
                            legend.position = c("topright", "right", "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "none"),
                            ...){

  if(class(x) != "decision_curve") stop("x must be an object of class 'decision_curve'-- plot_roc_components is only functional for one decision curve at a time.")

  legend.position <- match.arg(legend.position)

  #if(missing(curve.names)) curve.names  <- NA
  if(missing(confidence.intervals)) confidence.intervals <- NA

  prepData <- preparePlotData(x = x,
                              curve.names = NA,
                              confidence.intervals = confidence.intervals)

  predictors <- prepData$predictors
  dc.data <- prepData$dc.data
  confidence.intervals <- prepData$confidence.intervals
  rm(prepData)

  #set some defaults if needed
  if(missing(xlim)) xlim = range(dc.data$thresholds)

  if(length(col) ==1 ) col <- rep(col, 2)

  if(missing(lwd)) lwd = 2
  if(length(lwd) ==1) lwd <- rep(lwd, length(predictors))
  if(length(lwd) == length(predictors)) lwd = c(lwd, 1, 1)

  if(missing(ylab)) ylab <-  paste("Number high risk (out of ", population.size,  ")", sep = "")

  if(missing(ylim)) ylim = c(0, population.size*1.05)

  plot_generic(xx = dc.data,
               predictors = predictors,
               value = "prob.high.risk",
               plotNew = TRUE,
               standardize = FALSE,
               confidence.intervals,
               cost.benefit.axis = cost.benefit.axis,
               cost.benefits = cost.benefits,
               n.cost.benefits = n.cost.benefits,
               cost.benefit.xlab = cost.benefit.xlab,
               xlab = xlab, ylab = ylab,
               col = col,
               lty = lty, lwd = lwd,
               xlim = xlim, ylim = ylim,
               legend.position = "none",
               population.size = population.size,
               ...) #add my own legend

  plot_generic(xx = dc.data,
               predictors = predictors,
               value = "DP",
               plotNew = FALSE,
               standardize = FALSE,
               confidence.intervals,
               cost.benefit.axis = cost.benefit.axis,
               cost.benefits = cost.benefits,
               n.cost.benefits = n.cost.benefits,
               cost.benefit.xlab = cost.benefit.xlab,
               xlab = xlab, ylab = ylab,
               col = col,
               lty = 2, lwd = lwd,
               xlim = xlim, ylim = ylim,
               legend.position = legend.position,
               lty.fpr = 0,
               lty.tpr = 0,
               tpr.fpr.legend = FALSE,
               impact.legend = TRUE,
               population.size = population.size,
               ...)

}
