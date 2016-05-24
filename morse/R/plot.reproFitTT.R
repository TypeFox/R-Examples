#' Plotting method for \code{reproFitTT} objects
#'
#' @param x an object of class \code{reproFitTT}
#' @param xlab a title for the \eqn{x}-label
#' @param ylab a title for the \eqn{y}-label
#' @param main main title for the plot
#' @param fitcol color used for the fitted curve
#' @param fitlty line type for the fitted curve
#' @param fitlwd width of the fitted curve
#' @param ci if \code{TRUE}, draws the 95 \% credible limits of the fitted curve
#' @param cicol color for the 95 \% credible limits of the fitted curve
#' @param cilty line type for the 95 \% credible limits of the fitted curve
#' @param cilwd width of the 95 \% credible limits of the fitted curve
#' @param addlegend if \code{TRUE}, adds a default legend to the plot
#' @param log.scale if \code{TRUE}, displays \eqn{x}-axis in log-scale
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param \dots Further arguments to be passed to generic methods.
#' @note When \code{style = "ggplot"}, the function calls package
#' \code{\link[ggplot2]{ggplot}} and returns an object of class \code{ggplot}.
#' @note For an example, see the paragraph on \code{\link{reproFitTT}}.
#' 
#' @import ggplot2
#' @import grDevices
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot axis legend lines par points polygon
#' segments title
#'
#' @keywords plot
#' 
#' @export
plot.reproFitTT <- function(x,
                            xlab = "Concentration",
                            ylab = "Nb of offspring per ind.day",
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
  # plot the fitted curve estimated by reproFitTT
  # INPUTS
  # - x:  reproFitTT object
  # - xlab : label x
  # - ylab : label y
  # - main : main title
  # - fitcol : color fitted curve
  # - fitlty : type line fitted curve
  # - fitlwd : width line fitted curve
  # - ci : credible interval, boolean
  # - cicol : color ci
  # - cilty : type line ci
  # - cilwd : width line ci
  # - addlegend : boolean
  # - log.scale : x log option
  # - style : generic ou ggplot
  # OUTPUT:
  # - plot of fitted regression

  # Selection of datapoints that can be displayed given the type of scale
  sel <- if(log.scale) x$dataTT$conc > 0 else TRUE

  dataTT <- x$dataTT[sel, ]
  dataTT$resp <- dataTT$Nreprocumul / dataTT$Nindtime
  transf_data_conc <- optLogTransform(log.scale, dataTT$conc)

  # Concentration values used for display in linear scale
  display.conc <- (function() {
    x <- optLogTransform(log.scale, dataTT$conc)
    s <- seq(min(x),max(x), length = 100)
    if(log.scale) exp(s) else s
  })()

  # Possibly log transformed concentration values for display
  curv_conc <- optLogTransform(log.scale, display.conc)

  curv_resp <- reproEvalFit(x, display.conc)

  CI <- if (ci) { CI <- reproLlmCI(x, display.conc) } else NULL

  if (style == "generic") {
    reproFitPlotGeneric(dataTT$conc, transf_data_conc, dataTT$resp,
                        curv_conc, curv_resp,
                        CI,
                        xlab, ylab, fitcol, fitlty, fitlwd,
                        main, addlegend,
                        cicol, cilty, cilwd)
  }
  else if (style == "ggplot") {
    reproFitPlotGG(dataTT$conc, transf_data_conc, dataTT$resp,
                   curv_conc, curv_resp,
                   CI,
                   xlab, ylab, fitcol, fitlty, fitlwd,
                   main, addlegend,
                   cicol, cilty, cilwd)
  }
  else stop("Unknown style")
}


reproEvalFit <- function(fit, x) {
  # eval the fitted function on x
  # INPUT :
  # - fit: reproFitTT object
  # - x: vector of concentrations
  # OUTPUT :
  # - fNcumulpidtheo

  res.M <- summary(fit$mcmc)

  # unlog parameters
  d <- res.M$quantiles["d", "50%"]
  b <- 10^res.M$quantiles["log10b", "50%"]
  e <- 10^res.M$quantiles["log10e", "50%"]
  fNcumulpidtheo <- d / (1 + ( x / e)^b) # mean curve equation

  return(fNcumulpidtheo)
}

#' @importFrom stats quantile rgamma
reproLlmCI <- function(fit, x) {
  # create the parameters for credible interval for the log logistic model
  # INPUT:
  # - fit : object of class reproFitTT
  # - x : vector of concentrations values (x axis)
  # OUTPUT:
  # - ci : credible limit

  mctot <- do.call("rbind", fit$mcmc)
  k <- nrow(mctot)
  # parameters
  d2 <- mctot[, "d"]
  log10b2 <- mctot[, "log10b"]
  b2 <- 10^log10b2
  log10e2 <- mctot[, "log10e"]
  e2 <- 10^log10e2

  # quantiles
  qinf95 = NULL
  qsup95 = NULL

  # poisson
  if (fit$model.label == "P") {
    for (i in 1:length(x)) {
      theomean <- d2 / (1 + (x[i] / e2)^(b2)) # mean curve
      # IC 95%
      qinf95[i] <- quantile(theomean, probs = 0.025, na.rm = TRUE)
      qsup95[i] <- quantile(theomean, probs = 0.975, na.rm = TRUE)
    }
  }

  # gamma poisson
  else if (fit$model.label == "GP") {
    # parameters
    log10omega2 <- mctot[, "log10omega"]
    omega2 <- 10^(log10omega2)

    for (i in 1:length(x)) {
      theomean <- d2 / (1 + (x[i] / e2)^(b2)) # mean curve
      theo <- rgamma(n = k, shape = theomean / omega2, rate = 1 / omega2)
      # IC 95%
      qinf95[i] <- quantile(theo, probs = 0.025, na.rm = TRUE)
      qsup95[i] <- quantile(theo, probs = 0.975, na.rm = TRUE)
    }
  }
  # values for CI
  ci <- list(qinf95 = qinf95,
             qsup95 = qsup95)

  return(ci)
}

reproFitPlotGenericNoCI <- function(data_conc, transf_data_conc, data_resp,
                                    curv_conc, curv_resp,
                                    xlab, ylab, fitcol, fitlty, fitlwd,
                                    main, addlegend) {
  # plot the fitted curve estimated by reproFitTT
  # with generic style without credible interval

  plot(transf_data_conc, data_resp,
       xlab = xlab,
       ylab = ylab,
       main = main,
       xaxt = "n",
       yaxt = "n",
       ylim = c(0, max(data_resp) + 0.2),
       type = "n")
  # axis
  axis(side = 2, at = pretty(c(0, max(data_resp))))
  axis(side = 1, at = transf_data_conc,
       labels = data_conc)

  # fitted curve
  lines(curv_conc, curv_resp, col = fitcol,
        lty = fitlty, lwd = fitlwd, type = "l")

  # legend
  if (addlegend) {
    legend("bottomleft",
           lty = fitlty,
           lwd = fitlwd,
           col = fitcol,
           legend = "loglogistic",
           bty = "n")
  }
}

reproFitPlotGenericCI <- function(data_conc, transf_data_conc, data_resp,
                                  curv_conc, curv_resp,
                                  CI,
                                  xlab, ylab, fitcol, fitlty, fitlwd,
                                  main, addlegend,
                                  cicol, cilty, cilwd) {
  # plot the fitted curve estimated by reproFitTT
  # with generic style with credible interval

  plot(transf_data_conc, data_resp,
       xlab = xlab,
       ylab = ylab,
       main = main,
       xaxt = "n",
       yaxt = "n",
       ylim = c(0, max(c(data_resp, CI[["qsup95"]])) + 0.01),
       type = "n")

  # axis
  axis(side = 2, at = pretty(c(0, max(CI[["qsup95"]]))))
  axis(side = 1,
       at = transf_data_conc,
       labels = data_conc)

  # Plotting the theoretical curve
  # CI ribbon + lines
  polygon(c(curv_conc, rev(curv_conc)), c(CI[["qinf95"]], rev(CI[["qsup95"]])),
          col = cicol)
  lines(curv_conc, CI[["qsup95"]], type = "l", col = cicol, lty = cilty,
        lwd = cilwd)
  lines(curv_conc, CI[["qinf95"]], type = "l", col = cicol, lty = cilty,
        lwd = cilwd)

  # fitted curve
  lines(curv_conc, curv_resp, col = fitcol,
        lty = fitlty, lwd = fitlwd, type = "l")

  # legend
  if(addlegend)
  legend("bottomleft",
         lty = c(fitlty, cilty),
         lwd = c(fitlwd, cilwd),
         col = c(fitcol, cicol),
         legend = c("loglogistic", "95% Credible limits"),
         bty = "n")
}

reproFitPlotGeneric <- function(data_conc, transf_data_conc, data_resp,
                                curv_conc, curv_resp,
                                CI,
                                xlab, ylab, fitcol, fitlty, fitlwd,
                                main, addlegend,
                                cicol, cilty, cilwd) {

  if(!is.null(CI)) reproFitPlotGenericCI(data_conc, transf_data_conc,
                                         data_resp,
                                         curv_conc, curv_resp,
                                         CI,
                                         xlab, ylab, fitcol, fitlty, fitlwd,
                                         main, addlegend,
                                         cicol, cilty, cilwd)
  else {
    reproFitPlotGenericNoCI(data_conc, transf_data_conc, data_resp,
                            curv_conc, curv_resp,
                            xlab, ylab, fitcol, fitlty, fitlwd,
                            main, addlegend)
  }
}

reproFitPlotGGNoCI <- function(curv, valCols,
                               fitlty, fitlwd, xlab, ylab, main) {
  plt_4 <- ggplot(curv) +
    geom_line(aes(conc, resp), curv,
              linetype = fitlty, size = fitlwd, color = valCols$cols2) +
    ylim(0, max(curv$resp) + 1) +
    labs(x = xlab, y = ylab) +
    ggtitle(main) + theme_minimal()

  return(plt_4)
}

reproFitPlotGGCI <- function(curv, CI, cicol, cilty, cilwd,
                             valCols, fitlty, fitlwd, xlab, ylab, main) {
  # IC
  cri <- data.frame(conc = curv$conc,
                           qinf95 = CI[["qinf95"]],
                           qsup95 = CI[["qsup95"]],
                           CI = "Credible limits")

  plt_3 <- ggplot(cri) +
    geom_line(data = cri, aes(conc, qinf95, color = CI),
              linetype = cilty, size = cilwd) +
    geom_line(data = cri, aes(conc, qsup95, color = CI),
              linetype = cilty, size = cilwd) +
    geom_ribbon(data = cri, aes(x = conc, ymin = qinf95,
                                ymax = qsup95), fill = valCols$cols3,
                col = valCols$cols3,
                alpha = 0.4) +
    scale_color_manual(name = "", values = valCols$cols3) +
    theme_minimal()

  # plot IC
  # final plot
  plt_4 <- ggplot(cri) +
    geom_line(data = cri, aes(conc, qinf95),
              linetype = cilty, size = cilwd, color = valCols$cols3) +
    geom_line(data = cri, aes(conc, qsup95),
              linetype = cilty, size = cilwd, color = valCols$cols3) +
    geom_ribbon(data = cri, aes(x = conc, ymin = qinf95,
                                ymax = qsup95), fill = valCols$cols3,
                col = valCols$cols3,
                alpha = 0.4) +
    geom_line(aes(conc, resp), curv,
              linetype = fitlty, size = fitlwd, color = valCols$cols2) +
    ylim(0, max(CI[["qsup95"]]) + 0.2) +
    labs(x = xlab, y = ylab) +
    ggtitle(main) + theme_minimal()

  return(list(plt_3 = plt_3,
              plt_4 = plt_4))
}

reproFitPlotGG <- function(data_conc, transf_data_conc, data_resp,
                           curv_conc, curv_resp,
                           CI,
                           xlab, ylab, fitcol, fitlty, fitlwd,
                           main, addlegend,
                           cicol, cilty, cilwd) {

  if (Sys.getenv("RSTUDIO") == "") {
    dev.new() # create a new page plot
    # when not use RStudio
  }

  # dataframes points (data) and curve (curv)
  curv <- data.frame(conc = curv_conc, resp = curv_resp, Line = "loglogistic")

  # colors
  valCols <- fCols(curv, fitcol, cicol, "repro")

  # curve (to create the legend)
    plt_2 <- ggplot(curv) +
    geom_line(data = curv, aes(conc, resp, color = Line),
              linetype = fitlty, size = fitlwd) +
    scale_color_manual(name = "", values = valCols$cols2) +
    theme_minimal()

  plt_4 <-
    if (! is.null(CI)) {
      reproFitPlotGGCI(curv, CI, cicol, cilty, cilwd,
                       valCols, fitlty, fitlwd, xlab, ylab, main)$plt_4
    } else {
      reproFitPlotGGNoCI(curv, valCols, fitlty, fitlwd,
                         xlab, ylab, main)
    }

  if (addlegend) {
    # create legends
    mylegend_2 <- legendGgplotFit(plt_2) # mean line legend

    plt_5 <- plt_4 + scale_x_continuous(breaks = transf_data_conc,
                                        labels = data_conc)

    if (is.null(CI)) {
      grid.arrange(plt_5, arrangeGrob(mylegend_2, nrow = 6),
                   ncol = 2, widths = c(6, 2))
    }
    else {
      plt_3 <- reproFitPlotGGCI(curv, CI, cicol, cilty, cilwd,
                                valCols, fitlty, fitlwd, xlab, ylab, main)$plt_3
      mylegend_3 <- legendGgplotFit(plt_3)
      grid.arrange(plt_5, arrangeGrob(mylegend_2, mylegend_3,
                                      nrow = 6), ncol = 2,
                   widths = c(6, 2))
    }
  }
  else { # no legend
    plt_5 <- plt_4 + scale_x_continuous(breaks = transf_data_conc,
                                        labels = data_conc)
    return(plt_5)
  }
}

