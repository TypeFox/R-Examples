# @file Plots.R
#
# Copyright 2015 Observational Health Data Sciences and Informatics
#
# This file is part of EmpiricalCalibration
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Create a forest plot
#'
#' @description
#' \code{plotForest} creates a forest plot of effect size estimates.
#'
#' @details
#' Creates a forest plot of effect size estimates (ratios). Estimates that are significantly different
#' from 1 (alpha = 0.05) are marked in orange, others are marked in blue.
#'
#'
#' @param logRr      A numeric vector of effect estimates on the log scale
#' @param seLogRr    The standard error of the log of the effect estimates. Hint: often the standard
#'                   error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                   estimate>))/qnorm(0.025)
#' @param names      A vector containing the names of the drugs or outcomes
#' @param xLabel     The label on the x-axis: the name of the effect estimate
#'
#' @return
#' A Ggplot object. Use the \code{ggsave} function to save to file.
#' @param fileName   Name of the file where the plot should be saved, for example 'plot.png'. See the
#'                   function \code{ggsave} in the ggplot2 package for supported file formats.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' plotForest(negatives$logRr, negatives$seLogRr, negatives$drugName)
#'
#' @export
plotForest <- function(logRr, seLogRr, names, xLabel = "Relative risk", fileName = NULL) {
  breaks <- c(0.25, 0.5, 1, 2, 4, 6, 8, 10)
  theme <- ggplot2::element_text(colour = "#000000", size = 6)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 5, hjust = 1)
  col <- c(rgb(0, 0, 0.8, alpha = 1), rgb(0.8, 0.4, 0, alpha = 1))
  colFill <- c(rgb(0, 0, 1, alpha = 0.5), rgb(1, 0.4, 0, alpha = 0.5))
  data <- data.frame(DRUG_NAME = as.factor(names),
                     logRr = logRr,
                     logLb95Rr = logRr + qnorm(0.025) *
    seLogRr, logUb95Rr = logRr + qnorm(0.975) * seLogRr)
  data$significant <- data$logLb95Rr > 0 | data$logUb95Rr < 0
  data$DRUG_NAME <- factor(data$DRUG_NAME, levels = rev(levels(data$DRUG_NAME)))
  plot <- with(data, {
    ggplot2::ggplot(data,
                    ggplot2::aes(x = DRUG_NAME,
                                 y = exp(logRr),
                                 ymin = exp(logLb95Rr),
                                 ymax = exp(logUb95Rr),
                                 colour = significant,
                                 fill = significant),
                    environment = environment()) + ggplot2::geom_hline(yintercept = breaks,
                                                                       colour = "#AAAAAA",
                                                                       lty = 1,
                                                                       size = 0.2) + ggplot2::geom_hline(yintercept = 1, size = 0.5) + ggplot2::geom_pointrange(shape = 23) + ggplot2::scale_colour_manual(values = col) + ggplot2::scale_fill_manual(values = colFill) + ggplot2::coord_flip(ylim = c(0.25, 10)) + ggplot2::scale_y_continuous(xLabel, trans = "log10", breaks = breaks, labels = breaks) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA), panel.grid.major = ggplot2::element_line(colour = "#EEEEEE"), axis.ticks = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(), axis.text.y = themeRA, axis.text.x = theme, legend.key = ggplot2::element_blank(), strip.text.x = theme, strip.background = ggplot2::element_blank(), legend.position = "none")
  })
  if (!is.null(fileName))
    ggplot2::ggsave(fileName, plot, width = 5, height = 2.5 + length(logRr) * 0.8, dpi = 400)
  return(plot)
}

logRrtoSE <- function(logRr, p, null) {
  sapply(logRr, function(logRr) {
    precision <- 0.001
    if (calibrateP(null, logRr, precision, pValueOnly = TRUE) > p)
      return(0)
    L <- 0
    H <- 100
    while (H >= L) {
      M <- L + (H - L)/2
      if (calibrateP(null, logRr, M, pValueOnly = TRUE) - p > precision)
        H <- M else if (p - calibrateP(null, logRr, M, pValueOnly = TRUE) > precision)
        L <- M else return(M)
    }
    return(L - 1)
  })
}

#' Plot the effect of the calibration
#'
#' @description
#' \code{plotCalibrationEffect} creates a plot showing the effect of the calibration.
#'
#' @details
#' Creates a plot with the effect estimate on the x-axis and the standard error on the y-axis.
#' Negative controls are shown as blue dots, positive controls as yellow diamonds. The area below the
#' dashed line indicated estimates with p < 0.05. The orange area indicates estimates with calibrated
#' p < 0.05.
#'
#' @param logRrNegatives     A numeric vector of effect estimates of the negative controls on the log
#'                           scale.
#' @param seLogRrNegatives   The standard error of the log of the effect estimates of the negative
#'                           controls.
#' @param logRrPositives     A numeric vector of effect estimates of the positive controls on the log
#'                           scale.
#' @param seLogRrPositives   The standard error of the log of the effect estimates of the positive
#'                           controls.
#' @param null               An object representing the fitted null distribution as created by the
#'                           \code{fitNull} function.
#' @param xLabel             The label on the x-axis: the name of the effect estimate.
#' @param fileName           Name of the file where the plot should be saved, for example 'plot.png'.
#'                           See the function \code{ggsave} in the ggplot2 package for supported file
#'                           formats.
#'
#' @return
#' A Ggplot object. Use the \code{ggsave} function to save to file.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' positive <- sccs[sccs$groundTruth == 1, ]
#' plotCalibrationEffect(negatives$logRr, negatives$seLogRr, positive$logRr, positive$seLogRr)
#'
#' @export
plotCalibrationEffect <- function(logRrNegatives,
                                  seLogRrNegatives,
                                  logRrPositives,
                                  seLogRrPositives,
                                  null = NULL,
                                  xLabel = "Relative risk",
                                  fileName = NULL) {
  if (is.null(null))
    null <- fitNull(logRrNegatives, seLogRrNegatives)
  x <- exp(seq(log(0.25), log(10), by = 0.01))
  y <- logRrtoSE(log(x), 0.05, null)
  seTheoretical <- sapply(x, FUN = function(x) {
    abs(log(x))/qnorm(0.975)
  })
  breaks <- c(0.25, 0.5, 1, 2, 4, 6, 8, 10)
  theme <- ggplot2::element_text(colour = "#000000", size = 12)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 12, hjust = 1)
  plot <- ggplot2::ggplot(data.frame(x, y, seTheoretical),
                          ggplot2::aes(x = x, y = y),
                          environment = environment()) +
          ggplot2::geom_vline(xintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.5) +
          ggplot2::geom_vline(xintercept = 1, size = 1) +
          ggplot2::geom_area(fill = rgb(1, 0.5, 0, alpha = 0.5),
                             color = rgb(1, 0.5, 0),
                             size = 1,
                             alpha = 0.5) +
          ggplot2::geom_area(ggplot2::aes(y = seTheoretical),
                             fill = rgb(0, 0, 0),
                             colour = rgb(0, 0, 0, alpha = 0.1),
                             alpha = 0.1) +
          ggplot2::geom_line(ggplot2::aes(y = seTheoretical),
                             colour = rgb(0, 0, 0),
                             linetype = "dashed",
                             size = 1,
                             alpha = 0.5) +
          ggplot2::geom_point(shape = 21,
                              ggplot2::aes(x, y),
                              data = data.frame(x = exp(logRrNegatives), y = seLogRrNegatives),
                              size = 2,
                              fill = rgb(0, 0, 1, alpha = 0.5),
                              colour = rgb(0, 0, 0.8)) +
          ggplot2::geom_hline(yintercept = 0) +
          ggplot2::scale_x_continuous(xLabel,
                                      trans = "log10",
                                      limits = c(0.25, 10),
                                      breaks = breaks,
                                      labels = breaks) +
          ggplot2::scale_y_continuous("Standard Error", limits = c(0, 1.5)) +
          ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                         panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA),
                         panel.grid.major = ggplot2::element_blank(),
                         axis.ticks = ggplot2::element_blank(),
                         axis.text.y = themeRA,
                         axis.text.x = theme,
                         legend.key = ggplot2::element_blank(),
                         strip.text.x = theme,
                         strip.background = ggplot2::element_blank(),
                         legend.position = "none")
  if (!missing(logRrPositives)) {
    plot <- plot + ggplot2::geom_point(shape = 23,
                                       ggplot2::aes(x, y),
                                       data = data.frame(x = exp(logRrPositives),
                                                         y = seLogRrPositives),
                                       size = 4,
                                       fill = rgb(1, 1, 0),
                                       alpha = 0.8)
  }
  if (!is.null(fileName))
    ggplot2::ggsave(fileName, plot, width = 6, height = 4.5, dpi = 400)
  return(plot)
}

#' Create a calibration plot
#'
#' @description
#' \code{plotCalibration} creates a plot showing the calibration of our calibration procedure
#'
#' @details
#' Creates a calibration plot showing the number of effects with p < alpha for every level of alpha.
#' The empirical calibration is performed using a leave-one-out design: The p-value of an effect is
#' computed by fitting a null using all other negative controls. Ideally, the calibration line should
#' approximate the diagonal. The plot shows both theoretical (traditional) and empirically calibrated
#' p-values.
#'
#' @param logRr      A numeric vector of effect estimates on the log scale
#' @param seLogRr    The standard error of the log of the effect estimates. Hint: often the standard
#'                   error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                   estimate>))/qnorm(0.025)
#' @param useMcmc    Use MCMC to estimate the calibrated P-value?
#' @param fileName   Name of the file where the plot should be saved, for example 'plot.png'. See the
#'                   function \code{ggsave} in the ggplot2 package for supported file formats.
#'
#' @return
#' A Ggplot object. Use the \code{ggsave} function to save to file.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' plotCalibration(negatives$logRr, negatives$seLogRr)
#'
#' @export
plotCalibration <- function(logRr, seLogRr, useMcmc = FALSE, fileName = NULL) {
  if (any(is.infinite(seLogRr))) {
    warning("Estimate(s) with infinite standard error detected. Removing before fitting null distribution")
    logRr <- logRr[!is.infinite(seLogRr)]
    seLogRr <- seLogRr[!is.infinite(seLogRr)]
  }
  if (any(is.infinite(logRr))) {
    warning("Estimate(s) with infinite logRr detected. Removing before fitting null distribution")
    seLogRr <- seLogRr[!is.infinite(logRr)]
    logRr <- logRr[!is.infinite(logRr)]
  }
  if (any(is.na(seLogRr))) {
    warning("Estimate(s) with NA standard error detected. Removing before fitting null distribution")
    logRr <- logRr[!is.na(seLogRr)]
    seLogRr <- seLogRr[!is.na(seLogRr)]
  }
  if (any(is.na(logRr))) {
    warning("Estimate(s) with NA logRr detected. Removing before fitting null distribution")
    seLogRr <- seLogRr[!is.na(logRr)]
    logRr <- logRr[!is.na(logRr)]
  }

  data <- data.frame(logRr = logRr, SE = seLogRr)
  data$Z <- data$logRr/data$SE
  data$P <- 2 * pmin(pnorm(data$Z), 1 - pnorm(data$Z))  # 2-sided p-value
  data$Y <- sapply(data$P, function(x) {
    sum(data$P < x)/nrow(data)
  })

  data$calibratedP <- vector(length = nrow(data))
  for (i in 1:nrow(data)) {
    dataLeaveOneOut <- data[seq(1, nrow(data)) != i, ]
    if (useMcmc) {
      null <- fitMcmcNull(dataLeaveOneOut$logRr, dataLeaveOneOut$SE)
    } else {
      null <- fitNull(dataLeaveOneOut$logRr, dataLeaveOneOut$SE)
    }
    data$calibratedP[i] <- calibrateP(null, data$logRr[i], data$SE[i], pValueOnly = TRUE)
  }
  data <- data[!is.na(data$calibratedP), ]
  data$AdjustedY <- sapply(data$calibratedP, function(x) {
    sum(data$calibratedP < x)/nrow(data)
  })

  catData <- data.frame(x = c(data$P, data$calibratedP),
                        y = c(data$Y, data$AdjustedY),
                        label = factor(c(rep("Theoretical", times = nrow(data)),
                                         rep("Empirical", times = nrow(data)))))
  catData$label <- factor(catData$label, levels = c("Empirical", "Theoretical"))

  names(catData) <- c("x", "y", "P-value calculation")

  breaks <- c(0, 0.25, 0.5, 0.75, 1)
  theme <- ggplot2::element_text(colour = "#000000", size = 10)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 10, hjust = 1)
  plot <- with(catData, {
    ggplot2::ggplot(catData,
                    ggplot2::aes(x = x,
                                 y = y,
                                 colour = `P-value calculation`,
                                 linetype = `P-value calculation`),
                    environment = environment()) + ggplot2::geom_vline(xintercept = breaks,
                                                                       colour = "#AAAAAA",
                                                                       lty = 1,
                                                                       size = 0.3) + ggplot2::geom_vline(xintercept = 0.05, colour = "#888888", linetype = "dashed", size = 1) + ggplot2::geom_hline(yintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.3) + ggplot2::geom_abline(colour = "#AAAAAA", lty = 1, size = 0.3) + ggplot2::geom_step(direction = "hv", size = 1) + ggplot2::scale_colour_manual(values = c(rgb(0, 0, 0), rgb(0, 0, 0), rgb(0.5, 0.5, 0.5))) + ggplot2::scale_linetype_manual(values = c("solid", "twodash")) + ggplot2::scale_x_continuous("Alpha", limits = c(0, 1), breaks = c(breaks, 0.05), labels = c("", ".25", ".50", ".75", "1", ".05")) + ggplot2::scale_y_continuous("Fraction with p < alpha", limits = c(0, 1), breaks = breaks, labels = c("0", ".25", ".50", ".75", "1")) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA), panel.grid.major = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.text.y = themeRA, axis.text.x = theme, strip.text.x = theme, strip.background = ggplot2::element_blank(), legend.position = "right")
  })
  if (!is.null(fileName))
    ggplot2::ggsave(fileName, plot, width = 6, height = 4.5, dpi = 400)
  return(plot)
}

#' Plot true and observed values
#'
#' @description
#' Plot true and observed values, for example from a simulation study.
#'
#' @details
#' Creates a forest plot of effect size estimates (ratios). Estimates that are significantly different
#' from the true value (alpha = 0.05) are marked in orange, others are marked in blue.
#'
#' @param logRr       A numeric vector of effect estimates on the log scale.
#' @param seLogRr     The standard error of the log of the effect estimates. Hint: often the standard
#'                    error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                    estimate>))/qnorm(0.025).
#' @param trueLogRr   A vector of the true effect sizes.
#' @param xLabel      The label on the x-axis: the name of the effect estimate.
#' @param fileName    Name of the file where the plot should be saved, for example 'plot.png'. See the
#'                    function \code{ggsave} in the ggplot2 package for supported file formats.
#'
#' @return
#' A Ggplot object. Use the \code{ggsave} function to save to file.
#'
#' @examples
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
#' plotTrueAndObserved(data$logRr, data$seLogRr, data$trueLogRr)
#'
#' @export
plotTrueAndObserved <- function(logRr,
                                seLogRr,
                                trueLogRr,
                                xLabel = "Relative risk",
                                fileName = NULL) {
  breaks <- c(0.25, 0.5, 1, 2, 4, 6, 8, 10)
  theme <- ggplot2::element_text(colour = "#000000", size = 6)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 5, hjust = 1)
  col <- c(rgb(0, 0, 0.8, alpha = 1), rgb(0.8, 0.4, 0, alpha = 1))
  colFill <- c(rgb(0, 0, 1, alpha = 0.5), rgb(1, 0.4, 0, alpha = 0.5))
  data <- data.frame(logRr = logRr,
                     logLb95Rr = logRr + qnorm(0.025) * seLogRr,
                     logUb95Rr = logRr + qnorm(0.975) * seLogRr,
                     trueLogRr = trueLogRr,
                     trueRr = round(exp(trueLogRr), 2))
  data$significant <- data$logLb95Rr > data$trueLogRr | data$logUb95Rr < data$trueLogRr
  data <- data[order(data$trueLogRr, data$logRr), ]
  data$order <- 1:nrow(data)
  coverage <- aggregate(!significant ~ trueRr, data = data, mean)
  names(coverage)[2] <- "coverage"
  plot <- with(data, {
    ggplot2::ggplot(data,
                    ggplot2::aes(x = exp(logRr),
                                 y = order,
                                 xmin = exp(logLb95Rr),
                                 xmax = exp(logUb95Rr),
                                 colour = significant,
                                 fill = significant),
                    environment = environment()) + ggplot2::geom_vline(xintercept = breaks,
                                                                       colour = "#AAAAAA",
                                                                       lty = 1,
                                                                       size = 0.2) + ggplot2::geom_errorbarh(ggplot2::aes(x = trueRr, xmax = trueRr, xmin = trueRr), height = 1, color = rgb(0, 0, 0), size = 1) + ggplot2::geom_errorbarh(height = 0) + ggplot2::geom_point(shape = 21, size = 1.5) + ggplot2::scale_colour_manual(values = col) + ggplot2::scale_fill_manual(values = colFill) + ggplot2::coord_cartesian(xlim = c(0.25, 10)) + ggplot2::scale_x_continuous(xLabel, trans = "log10", breaks = breaks, labels = breaks) + ggplot2::facet_grid(trueRr ~ ., scales = "free_y", space = "free") + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA), panel.grid.major = ggplot2::element_line(colour = "#EEEEEE"), axis.ticks = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), axis.text.x = theme, legend.key = ggplot2::element_blank(), strip.text.y = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), legend.position = "none")
  })
  if (!is.null(fileName))
    ggplot2::ggsave(fileName, plot, width = 5, height = 7, dpi = 400)
  return(plot)
}

#' Plot the coverage
#'
#' @details
#' Plot the fractions of estimates where the true effect size is below, above or within the confidence
#' interval, for one or more true effect sizes.
#'
#' @param logRr       A numeric vector of effect estimates on the log scale
#' @param seLogRr     The standard error of the log of the effect estimates. Hint: often the standard
#'                    error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                    estimate>))/qnorm(0.025)
#' @param trueLogRr   A vector of the true effect sizes
#' @param region      Size of the confidence interval. Default is .95 (95 percent).
#' @param fileName    Name of the file where the plot should be saved, for example 'plot.png'. See the
#'                    function \code{ggsave} in the ggplot2 package for supported file formats.
#'
#' @examples
#' data <- simulateControls(n = 50 * 3, mean = 0, sd = 0.15, trueLogRr = log(c(1, 2, 4)))
#' plotCoverage(data$logRr, data$seLogRr, data$trueLogRr)
#'
#' @export
plotCoverage <- function(logRr, seLogRr, trueLogRr, region = 0.95, fileName = NULL) {
  data <- data.frame(logRr = logRr,
                     logLb95Rr = logRr + qnorm((1 - region)/2) * seLogRr,
                     logUb95Rr = logRr + qnorm(1 - (1 - region)/2) * seLogRr,
                     trueLogRr = trueLogRr,
                     trueRr = round(exp(trueLogRr), 2))
  vizD <- data.frame()
  for (trueRr in unique(data$trueRr)) {
    subset <- data[data$trueRr == trueRr, ]
    d <- data.frame(trueRr = trueRr, group = c("Below CI",
                                               "Within CI",
                                               "Above CI"), fraction = 0, pos = 0)
    d$fraction[1] <- mean(subset$trueLogRr < subset$logLb95Rr)
    d$fraction[2] <- mean(subset$trueLogRr >= subset$logLb95Rr & subset$trueLogRr <= subset$logUb95Rr)
    d$fraction[3] <- mean(subset$trueLogRr > subset$logUb95Rr)
    d$pos[1] <- d$fraction[1]/2
    d$pos[2] <- d$fraction[1] + (d$fraction[2]/2)
    d$pos[3] <- d$fraction[1] + d$fraction[2] + (d$fraction[3]/2)
    vizD <- rbind(vizD, d)
  }
  vizD$pos <- sapply(vizD$pos, function(x) {
    min(max(x, 0.05), 0.95)
  })

  vizD$label <- paste(round(100 * vizD$fraction), "%", sep = "")
  vizD$group <- factor(vizD$group, levels = c("Below CI", "Within CI", "Above CI"))
  theme <- ggplot2::element_text(colour = "#000000", size = 10)
  plot <- with(vizD, {
    ggplot2::ggplot(vizD, ggplot2::aes(x = as.factor(trueRr),
                                       y = fraction)) + ggplot2::geom_bar(ggplot2::aes(fill = group),
                                                                                                 stat = "identity",
                                                                                                 position = "stack",
                                                                                                 alpha = 0.8) + ggplot2::scale_fill_manual(values = c("#174a9f",
                                                                                                                                                      "#f9dd75",
                                                                                                                                                      "#f15222")) + ggplot2::geom_text(ggplot2::aes(label = label, y = pos), size = 3) + ggplot2::scale_x_discrete("True relative risk") + ggplot2::scale_y_continuous("Coverage") + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA), panel.grid.major = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), axis.text.x = theme, legend.key = ggplot2::element_blank(), legend.position = "right")
  })
  if (!is.null(fileName))
    ggplot2::ggsave(fileName, plot, width = 5, height = 3.5, dpi = 400)
  return(plot)
}


#' Plot the MCMC trace
#'
#' @details
#' Plot the trace of the MCMC for diagnostics purposes.
#'
#' @param mcmcNull   An object of type \code{mcmcNull} as generated using the \code{fitMcmcNull}
#'                   function.
#' @param fileName   Name of the file where the plot should be saved, for example 'plot.png'. See the
#'                   function \code{ggsave} in the ggplot2 package for supported file formats.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' null <- fitMcmcNull(negatives$logRr, negatives$seLogRr)
#' plotMcmcTrace(null)
#'
#' @export
plotMcmcTrace <- function(mcmcNull, fileName = NULL) {
  mcmc <- attr(mcmcNull, "mcmc")
  dataMean <- data.frame(x = 1:nrow(mcmc$chain), trace = as.numeric(ts(mcmc$chain[,
                                                                       1])), var = "Mean")
  dataPrecision <- data.frame(x = 1:nrow(mcmc$chain), trace = as.numeric(ts(mcmc$chain[,
                                                                            2])), var = "Precision")
  data <- rbind(dataMean, dataPrecision)
  plot <- with(data, {
    ggplot2::ggplot(data, ggplot2::aes(x = x,
                                       y = trace)) + ggplot2::geom_line(alpha = 0.7) + ggplot2::scale_x_continuous("Iterations") + ggplot2::facet_grid(var ~ ., scales = "free") + ggplot2::theme(axis.title.y = ggplot2::element_blank())
  })
  if (!is.null(fileName))
    ggplot2::ggsave(fileName, plot, width = 5, height = 3.5, dpi = 400)
  return(plot)
}
