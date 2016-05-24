#' Plotting method for \code{survFitTT} objects
#'
#' This function plots exposure-response fits for target time survival
#' analysis (a.k.a. \code{survFitTT} objects).
#'
#' The fitted curve represents the \strong{estimated survival rate} after
#' the target time has passed as a function of the concentration of pollutant;
#' the black dots depict the \strong{observed survival rate} at each tested
#' concentration. Note that since our model does not take inter-replicate
#' variability into consideration, replicates are systematically pooled in this
#' plot. When \code{ci = TRUE}, the function plots both credible intervals for
#' the estimated survival rate (by default the red area around the fitted
#' curve) and confidence intervals for the observed survival rate (as black
#' error bars). Both types of intervals are taken at the same level. Typically
#' a good fit is expected to display a large overlap between the two intervals.
#'
#' @param x an object of class \code{survFitTT}
#' @param xlab a title for the \eqn{x}-axis
#' @param ylab a title for the \eqn{y}-axis
#' @param main main title for the plot
#' @param fitcol color of the fitted curve
#' @param fitlty line type of the fitted curve
#' @param fitlwd width of the fitted curve
#' @param ci if \code{TRUE}, draws the 95 \% confidence interval on observed data
#' @param cicol color of the 95 \% confidence interval limits
#' @param cilty line type for the 95 \% confidence interval limits
#' @param cilwd width of the 95 \% confidence interval limits
#' @param addlegend if \code{TRUE}, adds a default legend to the plot
#' @param log.scale if \code{TRUE}, displays \eqn{x}-axis in log scale
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param \dots Further arguments to be passed to generic methods.
#' @note When \code{style = "ggplot"}, the function calls package
#' \code{\link[ggplot2]{ggplot}} and returns an object of class \code{ggplot}.
#' @note For an example, see the paragraph on \code{\link{reproFitTT}}.
#'
#' @keywords plot
#'
#' @import grDevices
#' @import ggplot2
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot axis legend lines par points polygon
#' @importFrom stats aggregate
#'
#' @export
plot.survFitTT <- function(x,
                           xlab = "Concentration",
                           ylab = "Survival rate",
                           main = NULL,
                           fitcol = "red",
                           fitlty = 1,
                           fitlwd = 1,
                           ci = FALSE,
                           cicol = "pink1",
                           cilty = 1,
                           cilwd = 1,
                           addlegend = FALSE,
                           log.scale = FALSE,
                           style = "generic", ...) {
  # plot the fitted curve estimated by survFitTT
  # INPUTS
  # - x:  survFitTt object
  # - xlab : label x
  # - ylab : label y
  # - main : main title
  # - fitcol : color fitted curve
  # - fitlty : type line fitted curve
  # - fitlwd : width line fitted curve
  # - ci : credible interval, boolean
  # - cicol : color ci ribbon
  # - cilty : type line ci ribbon
  # - cilwd : width line ci ribbon
  # - addlegend : boolean
  # - log.scale : x log option
  # - style : generic or ggplot
  # OUTPUT:
  # - plot of fitted regression

  # Selection of datapoints that can be displayed given the type of scale
  sel <- if(log.scale) x$dataTT$conc > 0 else TRUE

  dataTT <- x$dataTT[sel, ]
  dataTT$resp <- dataTT$Nsurv / dataTT$Ninit
  # data points are systematically pooled, since our model does not
  # take individual variation into account
  dataTT <- aggregate(resp ~ conc, dataTT, mean)
  transf_data_conc <- optLogTransform(log.scale, dataTT$conc)

  # Concentration values used for display in linear scale
  display.conc <- (function() {
    x <- optLogTransform(log.scale, dataTT$conc)
    s <- seq(min(x),max(x), length = 100)
    if(log.scale) exp(s) else s
  })()

  # Possibly log transformed concentration values for display
  curv_conc <- optLogTransform(log.scale, display.conc)

  curv_resp <- survEvalFit(x, display.conc)

  conf.int <- if(ci) { survLlbinomConfInt(x, log.scale) } else NULL
  cred.int <- if(ci) { survMeanCredInt(x, display.conc) } else NULL

  if (style == "generic") {
    survFitPlotGeneric(x,
                       dataTT$conc, transf_data_conc, dataTT$resp,
                       curv_conc, curv_resp,
                       conf.int, cred.int, log.scale,
                       xlab, ylab, fitcol, fitlty, fitlwd,
                       main, addlegend,
                       cicol, cilty, cilwd)
  }
  else if (style == "ggplot") {
    survFitPlotGG(x,
                  dataTT$conc, transf_data_conc, dataTT$resp,
                  curv_conc, curv_resp,
                  conf.int, cred.int,
                  xlab, ylab, fitcol, fitlty, fitlwd,
                  main, addlegend,
                  cicol, cilty, cilwd / 2)
  }
  else stop("Unknown style")
}

survEvalFit <- function(fit, x) {
  # eval the fitted function on x
  # INPUT :
  # - fit: survFitTT object
  # - x: vector of concentrations
  # OUTPUT :
  # - fNsurvtheo

  res.M <- summary(fit$mcmc)

  # unlog parameters
  b <- 10^res.M$quantiles["log10b", "50%"]
  e <- 10^res.M$quantiles["log10e", "50%"]

  if (fit$det.part == "loglogisticbinom_3") {
    d <- res.M$quantiles["d", "50%"]
    d / (1 + (x / e)^b) # mean curve equation 3 parameters
  } else {
    1 / (1 + (x / e)^b) # mean curve equation 2 parameters
  }
}

#' @importFrom stats aggregate binom.test
survLlbinomConfInt <- function(x, log.scale) {
  # create confidente interval on observed data for the log logistic
  # binomial model by a binomial test
  # INPUT:
  # - x : object of class survFitTT
  # - log.scale : boolean
  # OUTPUT:

  # - ci : confidente interval
  x <- cbind(aggregate(Nsurv ~ time + conc, x$dataTT, sum),
             Ninit = aggregate(Ninit ~ time + conc, x$dataTT, sum)$Ninit)

  ci <- apply(x, 1, function(x) {
    binom.test(x["Nsurv"], x["Ninit"])$conf.int
  })
  rownames(ci) <- c("qinf95", "qsup95")
  colnames(ci) <- x$conc

  if (log.scale) ci <- ci[ ,colnames(ci) != 0]

  return(ci)
}

#' @importFrom stats quantile
survMeanCredInt <- function(fit, x) {
  # create the parameters for credible interval for the log logistic binomial
  # model
  # INPUT:
  # - fit : object of class survFitTT
  # - x : vector of concentrations values (x axis)
  # OUTPUT:
  # - ci : credible limit

  mctot <- do.call("rbind", fit$mcmc)
  k <- nrow(mctot)
  # parameters
  if (fit$det.part == "loglogisticbinom_3") {
    d2 <- mctot[, "d"]
  }
  log10b2 <- mctot[, "log10b"]
  b2 <- 10^log10b2
  log10e2 <- mctot[, "log10e"]
  e2 <- 10^log10e2

  # quantiles
  qinf95 = NULL
  qsup95 = NULL

  for (i in 1:length(x)) {
    # llbinom 2 parameters
    if (fit$det.part == "loglogisticbinom_2") {
      theomean <- 1 / (1 + (x[i] / e2)^(b2)) # mean curve
    }

    # llbinom 3 parameters
    else if (fit$det.part == "loglogisticbinom_3") {
      theomean <- d2 / (1 + (x[i] / e2)^(b2)) # mean curve
    }
    # IC 95%
    qinf95[i] <- quantile(theomean, probs = 0.025, na.rm = TRUE)
    qsup95[i] <- quantile(theomean, probs = 0.975, na.rm = TRUE)
  }

  # values for CI
  ci <- list(qinf95 = qinf95,
             qsup95 = qsup95)

  return(ci)
}

survFitPlotGenericNoCredInt <- function(x,
                                   data_conc, transf_data_conc, data_resp,
                                   curv_conc, curv_resp,
                                   xlab, ylab, fitcol, fitlty, fitlwd,
                                   main, addlegend)
{
  # plot the fitted curve estimated by survFitTT
  # with generic style without credible interval

  plot(transf_data_conc, data_resp,
       xlab = xlab,
       ylab = ylab,
       main = main,
       pch  = 16,
       xaxt = "n",
       yaxt = "n",
       ylim = c(0, 1.05))

  # axis
  axis(side = 2, at = pretty(c(0, 1)))
  axis(side = 1,
       at = transf_data_conc,
       labels = data_conc)

  # fitted curve
  lines(curv_conc, curv_resp, col = fitcol,
        lty = fitlty, lwd = fitlwd, type = "l")

  # legend
  if (addlegend) {
    legend("bottomleft",
           pch = c(16, NA),
           lty = c(NA, fitlty),
           lwd = c(NA, fitlwd),
           col = c(1, fitcol),
           legend = c("Observed values", x$det.part),
           bty = "n")
  }
}

survFitPlotGenericCredInt <- function(x,
                                 data_conc, transf_data_conc, data_resp,
                                 curv_conc, curv_resp,
                                 conf.int, cred.int, log.scale,
                                 xlab, ylab, fitcol, fitlty, fitlwd,
                                 main, addlegend,
                                 cicol, cilty, cilwd)
{
  # plot the fitted curve estimated by survFitTT
  # with generic style with credible interval
  plot(transf_data_conc, data_resp,
       xlab = xlab,
       ylab = ylab,
       main = main,
       xaxt = "n",
       yaxt = "n",
       ylim = c(0, max(conf.int["qsup95",]) + 0.01),
       type = "n")

  # axis
  axis(side = 2, at = pretty(c(0, max(conf.int["qsup95",]))))
  axis(side = 1,
       at = transf_data_conc,
       labels = data_conc)

  # Plotting the theoretical curve
  # CI ribbon + lines
  polygon(c(curv_conc, rev(curv_conc)), c(cred.int[["qinf95"]], rev(cred.int[["qsup95"]])),
          col = cicol)
  lines(curv_conc, cred.int[["qsup95"]], type = "l", col = cicol, lty = cilty,
        lwd = cilwd)
  lines(curv_conc, cred.int[["qinf95"]], type = "l", col = cicol, lty = cilty,
        lwd = cilwd)

  # segment CI

  segments(transf_data_conc, data_resp,
           transf_data_conc, conf.int["qsup95", ])

  Bond <- if (log.scale) {
    0.03 * (max(transf_data_conc) - min(transf_data_conc))
  } else {
    0.03 * (max(transf_data_conc) - min(transf_data_conc[which(transf_data_conc != 0)]))
  }

  segments(transf_data_conc - Bond,
           conf.int["qsup95", ],
           transf_data_conc + Bond,
           conf.int["qsup95", ])

  segments(transf_data_conc, data_resp,
           transf_data_conc, conf.int["qinf95", ])

  segments(transf_data_conc - Bond,
           conf.int["qinf95", ],
           transf_data_conc + Bond,
           conf.int["qinf95", ])

  # points
  points(transf_data_conc, data_resp, pch = 16)
  # fitted curve
  lines(curv_conc, curv_resp, col = fitcol,
        lty = fitlty, lwd = fitlwd, type = "l")

  # legend
  if (addlegend) {
    legend("bottomleft", pch = c(16, NA, NA, NA),
           lty = c(NA, 1, cilty, fitlty),
           lwd = c(NA, 1, cilwd, fitlwd),
           col = c(1, 1, cicol, fitcol),
           legend = c("Observed values", "Confidence interval",
                      "Credible limits", x$det.part),
           bty = "n")
  }
}


survFitPlotGeneric <- function(x,
                               data_conc, transf_data_conc, data_resp,
                               curv_conc, curv_resp,
                               conf.int, cred.int, log.scale,
                               xlab, ylab, fitcol, fitlty, fitlwd,
                               main, addlegend,
                               cicol, cilty, cilwd) {


  if(!is.null(conf.int)) survFitPlotGenericCredInt(x,
                                             data_conc, transf_data_conc, data_resp,
                                             curv_conc, curv_resp,
                                             conf.int, cred.int, log.scale,
                                             xlab, ylab, fitcol, fitlty, fitlwd,
                                             main, addlegend,
                                             cicol, cilty, cilwd)
  else {
    survFitPlotGenericNoCredInt(x,
                                data_conc, transf_data_conc, data_resp,
                                curv_conc, curv_resp,
                                xlab, ylab, fitcol, fitlty, fitlwd,
                                main, addlegend)
  }
}


survFitPlotGGNoCredInt <- function(data, curv, valCols,
                                   fitlty, fitlwd, xlab, ylab, main) {
  plt_4 <- ggplot(data) +
    geom_point(data = data, aes(transf_conc, resp)) +
    geom_line(aes(conc, resp), curv,
              linetype = fitlty, size = fitlwd, color = valCols$cols2) +
    scale_color_discrete(guide = "none") +
    ylim(0, 1) +
    labs(x = xlab, y = ylab) +
    ggtitle(main) + theme_minimal()

  return(plt_4)
}

#' @importFrom grid arrow unit
survFitPlotGGCredInt <- function(x, data, curv, conf.int, cred.int, cilty, cilwd,
                                 valCols, fitlty, fitlwd, xlab, ylab, main) {
  # IC
  data.three <- data.frame(conc = data$transf_conc,
                           qinf95 = conf.int["qinf95",],
                           qsup95 = conf.int["qsup95",],
                           Conf.Int = "Confidence interval")
  data.four <- data.frame(conc = curv$conc,
                          qinf95 = cred.int[["qinf95"]],
                          qsup95 = cred.int[["qsup95"]],
                          Cred.Lim = "Credible limits")

  plt_3 <- ggplot(data) +
    geom_segment(aes(x = conc, xend = conc, y = qinf95, yend = qsup95,
                     linetype = Conf.Int),
                 arrow = arrow(length = unit(0.25 , "cm"), angle = 90,
                              ends = "both"), data.three,
                 color = valCols$cols3) +
    scale_linetype(name = "") +
    theme_minimal()

  plt_32 <- ggplot(data) +
    geom_line(data = data.four, aes(conc, qinf95, color = Cred.Lim),
              linetype = cilty, size = cilwd) +
    geom_line(data = data.four, aes(conc, qsup95, color = Cred.Lim),
              linetype = cilty, size = cilwd) +
    scale_color_manual(name = "", values = valCols$cols4) +
    geom_ribbon(data = data.four, aes(x = conc, ymin = qinf95,
                                      ymax = qsup95),
                fill = valCols$cols4, col = valCols$cols4, alpha = 0.4) +
    theme_minimal()

  # plot IC
  # final plot
  plt_4 <- ggplot(data) +
    geom_segment(aes(x = conc, xend = conc, y = qinf95, yend = qsup95,
                     color = Conf.Int),
                 arrow = arrow(length = unit(0.25 , "cm"), angle = 90,
                               ends = "both"),
                 data.three, linetype = cilty,
                 size = cilwd, col = valCols$cols3) +
    geom_line(data = data.four, aes(conc, qinf95, color = Cred.Lim),
              linetype = cilty, size = cilwd) +
    geom_line(data = data.four, aes(conc, qsup95, color = Cred.Lim),
              linetype = cilty, size = cilwd) +
    geom_ribbon(data = data.four, aes(x = conc, ymin = qinf95,
                                      ymax = qsup95),
                fill = valCols$cols4,
                col = valCols$cols4, alpha = 0.4) +
    geom_point(data = data, aes(transf_conc, resp)) +
    geom_line(aes(conc, resp), curv, linetype = fitlty,
              size = fitlwd, color = valCols$cols2) +
    scale_color_discrete(guide = "none") +
    ylim(0, 1) +
    labs(x = xlab, y = ylab) +
    ggtitle(main) + theme_minimal()

  return(list(plt_3 = plt_3,
              plt_32 = plt_32,
              plt_4 = plt_4))
}

survFitPlotGG <- function(x,
                          data_conc, transf_data_conc, data_resp,
                          curv_conc, curv_resp,
                          conf.int, cred.int,
                          xlab, ylab, fitcol, fitlty, fitlwd,
                          main, addlegend,
                          cicol, cilty, cilwd) {


  if (Sys.getenv("RSTUDIO") == "") {
    dev.new() # create a new page plot
    # when not use RStudio
  }

  # dataframes points (data) and curve (curv)
  data <- data.frame(conc = data_conc, transf_conc = transf_data_conc,
                     resp = data_resp, Points = "Observed values")
  curv <- data.frame(conc = curv_conc, resp = curv_resp, Line = "loglogistic")

  # colors
  valCols <- fCols(data, fitcol, cicol, "surv")

  # points (to create the legend)
  plt_1 <- ggplot(data) +
    geom_point(data = data, aes(transf_conc, resp, color = Points)) +
    scale_color_manual(name = "", values = valCols$cols1) +
    theme_minimal()

  # curve (to create the legend)
  plt_2 <- ggplot(data) +
    geom_line(data = curv, aes(conc, resp, color = Line),
              linetype = fitlty, size = fitlwd) +
    scale_color_manual(name = "", values = valCols$cols2) +
    theme_minimal()

  plt_4 <-
    if (! is.null(conf.int)) {
      survFitPlotGGCredInt(x, data, curv, conf.int, cred.int, cilty, cilwd,
                      valCols, fitlty, fitlwd, xlab, ylab, main)$plt_4
    } else {
      survFitPlotGGNoCredInt(data, curv, valCols, fitlty, fitlwd,
                        xlab, ylab, main)
    }

  if (addlegend) { # legend yes
    # create legends
    mylegend_1 <- legendGgplotFit(plt_1) # points legend
    mylegend_2 <- legendGgplotFit(plt_2) # mean line legend

    plt_5 <- plt_4 + scale_x_continuous(breaks = data$transf_conc,
                                        labels = data$conc)

    if (is.null(conf.int)) {
      grid.arrange(plt_5, arrangeGrob(mylegend_1, mylegend_2, nrow = 6),
                   ncol = 2, widths = c(6, 2))
    }
    else {
      plt_3 <- survFitPlotGGCredInt(x, data, curv, conf.int, cred.int, cilty, cilwd,
                               valCols, fitlty, fitlwd, xlab, ylab, main)$plt_3
      plt_32 <- survFitPlotGGCredInt(x, data, curv, conf.int, cred.int, cilty, cilwd,
                                valCols, fitlty, fitlwd, xlab, ylab, main)$plt_32
      mylegend_3 <- legendGgplotFit(plt_3)
      mylegend_32 <- legendGgplotFit(plt_32)
      grid.arrange(plt_5, arrangeGrob(mylegend_1, mylegend_3, mylegend_32,
                                      mylegend_2, nrow = 6), ncol = 2,
                   widths = c(6, 2))
    }
  }
  else { # no legend
    plt_5 <- plt_4 + scale_x_continuous(breaks = data$transf_conc,
                                        labels = data$conc)
    return(plt_5)
  }
}

