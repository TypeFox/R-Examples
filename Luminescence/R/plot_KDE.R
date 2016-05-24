#' Plot kernel density estimate with statistics
#'
#' Plot a kernel density estimate of measurement values in combination with the
#' actual values and associated error bars in ascending order. Optionally,
#' statistical measures such as mean, median, standard deviation, standard
#' error and quartile range can be provided visually and numerically.
#'
#' The function allow passing several plot arguments, such as \code{main},
#' \code{xlab}, \code{cex}. However, as the figure is an overlay of two
#' separate plots, \code{ylim} must be specified in the order: c(ymin_axis1,
#' ymax_axis1, ymin_axis2, ymax_axis2) when using the cumulative values plot
#' option. Similarly, if other than the default colours are desired, the
#' argument col must be provided with colours in the following order:
#' probability density function, De values, De error bars, sd or qr polygon.
#' The line type (\code{lty}) for additional measures of centrality will cycle
#' through the default values (1, 2, ...) by default, i.e. KDE line solid,
#' further vertical lines dashed, dotted, dash-dotted and so on. To change this
#' behaviour specify the desired order of line types (e.g. \code{lty = c(1, 3,
#' 2, 5)}). See examples for some further explanations. For details on the
#' calculation of the bin-width (parameter \code{bw}) see
#' \code{\link{density}}.\cr\cr 
#' A statistic summary, i.e. a collection of statistic measures of 
#' centrality and dispersion (and further measures) can be added by specifying 
#' one or more of the following keywords: \code{"n"} (number of samples),
#' \code{"mean"} (mean De value), \code{"mean.weighted"} (error-weighted mean),
#' \code{"median"} (median of the De values), \code{"sdrel"} (relative standard
#' deviation in percent), \code{"sdrel.weighted"} (error-weighted relative 
#' standard deviation in percent), \code{"sdabs"} (absolute standard deviation),
#' \code{"sdabs.weighted"} (error-weighted absolute standard deviation), 
#' \code{"serel"} (relative standard error), \code{"serel.weighted"} (
#' error-weighted relative standard error), \code{"seabs"} (absolute standard
#' error), \code{"seabs.weighted"} (error-weighted absolute standard error), 
#' \code{"in.2s"} (percent of samples in 2-sigma range),
#' \code{"kurtosis"} (kurtosis) and \code{"skewness"} (skewness).
#'
#' @param data \code{\link{data.frame}} or \code{\linkS4class{RLum.Results}}
#' object (required): for \code{data.frame}: two columns: De
#' (\code{values[,1]}) and De error (\code{values[,2]}). For plotting multiple
#' data sets, these must be provided as \code{list} (e.g. \code{list(dataset1,
#' dataset2)}).
#' @param na.rm \code{\link{logical}} (with default): exclude NA values
#' from the data set prior to any further operations.
#' @param weights \code{\link{logical}} (with default): calculate the KDE with
#' De-errors as weights. Attention, using errors as weights will result in a
#' plot similar to a a probability density plot, with all ambiguities related
#' to this plot type!
#' @param values.cumulative \code{\link{logical}} (with default): show
#' cumulative individual data.
#' @param centrality \code{\link{character}}: measure(s) of centrality, used
#' for plotting vertical lines of the respective measure. Can be one out of
#' \code{"mean"}, \code{"median"}, \code{"mean.weighted"},
#' \code{"median.weighted"} and \code{"kdemax"}.
#' @param dispersion \code{\link{character}}: measure of dispersion, used for
#' drawing the polygon that depicts the dose distribution. One out of
#' \code{"sd"} (standard deviation),\code{"2sd"} (2 standard deviations)
#' \code{"qr"} (quartile range).
#' @param summary \code{\link{character}} (optional): add statistic measures of 
#' centrality and dispersion to the plot. Can be one or more of several 
#' keywords. See details for available keywords.
#' @param summary.pos \code{\link{numeric}} or \code{\link{character}} (with
#' default): optional position coordinates or keyword (e.g. \code{"topright"})
#' for the statistical summary. Alternatively, the keyword \code{"sub"} may be
#' specified to place the summary below the plot header. However, this latter
#' option in only possible if \code{mtext} is not used. In case of coordinate
#' specification, y-coordinate refers to the right y-axis.
#' @param polygon.col \code{\link{character}} or \code{\link{numeric}} (with
#' default): colour of the polygon showing the dose dispersion around the
#' central value. Only relevant if \code{dispersion} is specified.
#' @param order \code{\link{logical}}: Order data in ascending order.
#' @param bw \code{\link{character}} (with default): bin-width, chose a numeric
#' value for manual setting.
#' @param output \code{\link{logical}}: Optional output of numerical plot
#' parameters. These can be useful to reproduce similar plots. Default is
#' \code{FALSE}.
#' @param \dots further arguments and graphical parameters passed to
#' \code{\link{plot}}.
#' @note The plot output is no 'PD' plot (cf. the discussion of Berger and
#' Galbraith in Ancient TL; see references)!
#' @section Function version: 3.5
#' @author Michael Dietze, GFZ Potsdam (Germany),\cr Sebastian Kreutzer,
#' IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' @seealso \code{\link{density}}, \code{\link{plot}}
#' @examples
#'
#' ## read example data set
#' data(ExampleData.DeValues, envir = environment())
#' ExampleData.DeValues <-
#'   Second2Gray(ExampleData.DeValues$BT998, c(0.0438,0.0019))
#'
#' ## create plot straightforward
#' plot_KDE(data = ExampleData.DeValues)
#'
#' ## create plot with logarithmic x-axis
#' plot_KDE(data = ExampleData.DeValues,
#'          log = "x")
#'
#' ## create plot with user-defined labels and axes limits
#' plot_KDE(data = ExampleData.DeValues,
#'          main = "Dose distribution",
#'          xlab = "Dose (s)",
#'          ylab = c("KDE estimate", "Cumulative dose value"),
#'          xlim = c(100, 250),
#'          ylim = c(0, 0.08, 0, 30))
#'
#' ## create plot with centrality lines and distribution polygons
#' plot_KDE(data = ExampleData.DeValues,
#'          ylim = c(0, 0.08, 0, 35),
#'          centrality = c("median", "mean"),
#'          dispersion = "sd",
#'          polygon.col = "lightblue")
#'
#' ## create plot with statistical summary below header
#' plot_KDE(data = ExampleData.DeValues,
#'          summary = c("n", "median", "skewness", "qr"))
#'
#' ## create plot with statistical summary as legend
#' plot_KDE(data = ExampleData.DeValues,
#'          summary = c("n", "mean", "sdrel", "seabs"),
#'          summary.pos = "topleft")
#'
#' ## split data set into sub-groups, one is manipulated, and merge again
#' data.1 <- ExampleData.DeValues[1:15,]
#' data.2 <- ExampleData.DeValues[16:25,] * 1.3
#' data.3 <- list(data.1, data.2)
#'
#' ## create plot with two subsets straightforward
#' plot_KDE(data = data.3)
#'
#' ## create plot with two subsets and summary legend at user coordinates
#' plot_KDE(data = data.3,
#'          summary = c("n", "median", "skewness"),
#'          summary.pos = c(110, 0.07),
#'          col = c("blue", "orange"))
#'
#' ## example of how to use the numerical output of the function
#' ## return plot output to draw a thicker KDE line
#' KDE <- plot_KDE(data = ExampleData.DeValues,
#'                 output = TRUE)
#'
#' ## read out coordinates of KDE graph
#' KDE.x <- KDE$De.density[[1]]$x
#' KDE.y <- KDE$De.density[[1]]$y
#'
#' ## transform y-values to right y-axis dimensions
#' KDE.y <- KDE.y / max(KDE.y) * (nrow(ExampleData.DeValues) - 1) + 1
#'
#' ## draw the KDE line
#' lines(x = KDE.x,
#'       y = KDE.y,
#'       lwd = 3)
#'
#' @export
plot_KDE <- function(
  data,
  na.rm = TRUE,
  weights = FALSE,
  values.cumulative = TRUE,
  centrality,
  dispersion,
  summary,
  summary.pos,
  polygon.col,
  order = TRUE,
  bw = "nrd0",
  output = FALSE,
  ...
) {

  ## check data and parameter consistency -------------------------------------

  ## Homogenise input data format
  if(is(data, "list") == FALSE) {data <- list(data)}

  ## check/adjust input data structure
  for(i in 1:length(data)) {

    if(is(data[[i]], "RLum.Results") == FALSE &
         is(data[[i]], "data.frame") == FALSE &
         is.numeric(data[[i]]) == FALSE) {
      stop(paste("[plot_KDE()] Input data format is neither",
                 "'data.frame', 'RLum.Results' nor 'numeric'"))
    } else {

      if(is(data[[i]], "RLum.Results") == TRUE) {
        data[[i]] <- get_RLum(data[[i]], "data")
      }

      if(length(data[[i]]) < 2) {
        data[[i]] <- cbind(data[[i]], rep(NA, length(data[[i]])))
      }
    }
  }

  ## check/set function parameters
  if(missing(summary) == TRUE) {
    summary <- ""
  }

  if(missing(summary.pos) == TRUE) {
    summary.pos <- "sub"
  }

  mtext <- ""

  if(missing(polygon.col) == TRUE) {
    polygon.col <- rep("grey80",
                       length(data))
  }

  if(missing(centrality) == TRUE) {
    centrality <- character(0)
  }

  if(missing(dispersion) == TRUE) {
    dispersion <- ""
  }

  ## data preparation steps ---------------------------------------------------

  ## optionally, count and exclude NA values and print result
  if(na.rm == TRUE) {
    for(i in 1:length(data)) {
      n.NA <- sum(is.na(data[[i]][,1]))
      if(n.NA == 1) {
        print(paste("1 NA value excluded from data set", i, "."))
      } else if(n.NA > 1) {
        print(paste(n.NA, "NA values excluded from data set", i, "."))
      }
      data[[i]] <- na.exclude(data[[i]])
    }
  }

  ## optionally, order data set
  if(order == TRUE) {
    for(i in 1:length(data)) {
      data[[i]] <- data[[i]][order(data[[i]][,1]),]
    }
  }

  ## create output variables
  ## calculate and paste statistical summary
  De.stats <- matrix(nrow = length(data), ncol = 18)
  colnames(De.stats) <- c("n",
                          "mean",
                          "mean.weighted",
                          "median",
                          "median.weighted",
                          "kde.max",
                          "sd.abs",
                          "sd.rel",
                          "se.abs",
                          "se.rel",
                          "q25",
                          "q75",
                          "skewness",
                          "kurtosis",
                          "sd.abs.weighted",
                          "sd.rel.weighted",
                          "se.abs.weighted",
                          "se.rel.weighted")
  De.density <- list(NA)

  ## loop through all data sets
  for(i in 1:length(data)) {
    statistics <- calc_Statistics(data[[i]])
    De.stats[i,1] <- statistics$weighted$n
    De.stats[i,2] <- statistics$unweighted$mean
    De.stats[i,3] <- statistics$weighted$mean
    De.stats[i,4] <- statistics$unweighted$median
    De.stats[i,5] <- statistics$unweighted$median
    De.stats[i,7] <- statistics$unweighted$sd.abs
    De.stats[i,8] <- statistics$unweighted$sd.rel
    De.stats[i,9] <- statistics$unweighted$se.abs
    De.stats[i,10] <- statistics$weighted$se.rel
    De.stats[i,11] <- quantile(data[[i]][,1], 0.25)
    De.stats[i,12] <- quantile(data[[i]][,1], 0.75)
    De.stats[i,13] <- statistics$unweighted$skewness
    De.stats[i,14] <- statistics$unweighted$kurtosis
    De.stats[i,15] <- statistics$weighted$sd.abs
    De.stats[i,16] <- statistics$weighted$sd.rel
    De.stats[i,17] <- statistics$weighted$se.abs
    De.stats[i,18] <- statistics$weighted$se.rel

    De.density[[length(De.density) + 1]] <- if(weights == TRUE) {
      density(data[[i]][,1],
              kernel = "gaussian",
              bw = bw,
              weights = data[[i]][,2] / sum(data[[i]][,2]))
    } else {
      density(data[[i]][,1],
              kernel = "gaussian",
              bw = bw)
    }
  }

  ## remove dummy list element
  De.density[[1]] <- NULL

  ## create global data set
  De.global <- data[[1]][,1]
  De.error.global <- data[[1]][,2]
  De.density.range <- matrix(nrow = length(data), ncol = 4)
  for(i in 1:length(data)) {
    ##global De and De.error vector
    De.global <- c(De.global, data[[i]][,1])
    De.error.global <- c(De.error.global, data[[i]][,2])

    ## density ranges
    De.density.range[i,1] <- min(De.density[[i]]$x)
    De.density.range[i,2] <- max(De.density[[i]]$x)
    De.density.range[i,3] <- min(De.density[[i]]$y)
    De.density.range[i,4] <- max(De.density[[i]]$y)

    ## position of maximum KDE value
    De.stats[i,6] <- De.density[[i]]$x[which.max(De.density[[i]]$y)]
  }

  ## Get global range of densities
  De.density.range <- c(min(De.density.range[,1]),
                        max(De.density.range[,2]),
                        min(De.density.range[,3]),
                        max(De.density.range[,4]))

  label.text = list(NA)

  if(summary.pos[1] != "sub") {
    n.rows <- length(summary)

    for(i in 1:length(data)) {
      stops <- paste(rep("\n", (i - 1) * n.rows), collapse = "")

      summary.text <- character(0)

      for(j in 1:length(summary)) {
        summary.text <- c(summary.text,
                          paste(
                            "",
                            ifelse("n" %in% summary[j] == TRUE,
                                   paste("n = ",
                                         De.stats[i,1],
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("mean" %in% summary[j] == TRUE,
                                   paste("mean = ",
                                         round(De.stats[i,2], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("mean.weighted" %in% summary[j] == TRUE,
                                   paste("weighted mean = ",
                                         round(De.stats[i,3], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("median" %in% summary[j] == TRUE,
                                   paste("median = ",
                                         round(De.stats[i,4], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("median.weighted" %in% summary[j] == TRUE,
                                   paste("weighted median = ",
                                         round(De.stats[i,5], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("kdemax" %in% summary[j] == TRUE,
                                   paste("kdemax = ",
                                         round(De.stats[i,6], 2),
                                         " \n ",
                                         sep = ""),
                                   ""),
                            ifelse("sdabs" %in% summary[j] == TRUE,
                                   paste("sd = ",
                                         round(De.stats[i,7], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("sdrel" %in% summary[j] == TRUE,
                                   paste("rel. sd = ",
                                         round(De.stats[i,8], 2), " %",
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("seabs" %in% summary[j] == TRUE,
                                   paste("se = ",
                                         round(De.stats[i,9], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("serel" %in% summary[j] == TRUE,
                                   paste("rel. se = ",
                                         round(De.stats[i,10], 2), " %",
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("skewness" %in% summary[j] == TRUE,
                                   paste("skewness = ",
                                         round(De.stats[i,13], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("kurtosis" %in% summary[j] == TRUE,
                                   paste("kurtosis = ",
                                         round(De.stats[i,14], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("in.2s" %in% summary[j] == TRUE,
                                   paste("in 2 sigma = ",
                                         round(sum(data[[i]][,7] > -2 &
                                                     data[[i]][,7] < 2) /
                                                 nrow(data[[i]]) * 100 , 1),
                                         " %",
                                         sep = ""),
                                   ""),
                            ifelse("sdabs.weighted" %in% summary[j] == TRUE,
                                   paste("abs. weighted sd = ",
                                         round(De.stats[i,15], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("sdrel.weighted" %in% summary[j] == TRUE,
                                   paste("rel. weighted sd = ",
                                         round(De.stats[i,16], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("seabs.weighted" %in% summary[j] == TRUE,
                                   paste("abs. weighted se = ",
                                         round(De.stats[i,17], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("serel.weighted" %in% summary[j] == TRUE,
                                   paste("rel. weighted se = ",
                                         round(De.stats[i,18], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            sep = ""))
      }

      summary.text <- paste(summary.text, collapse = "")

      label.text[[length(label.text) + 1]] <- paste(stops,
                                                    summary.text,
                                                    stops,
                                                    sep = "")
    }
  } else {
    for(i in 1:length(data)) {

      summary.text <- character(0)

      for(j in 1:length(summary)) {
        summary.text <- c(summary.text,
                          ifelse("n" %in% summary[j] == TRUE,
                                 paste("n = ",
                                       De.stats[i,1],
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("mean" %in% summary[j] == TRUE,
                                 paste("mean = ",
                                       round(De.stats[i,2], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("mean.weighted" %in% summary[j] == TRUE,
                                 paste("weighted mean = ",
                                       round(De.stats[i,3], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("median" %in% summary[j] == TRUE,
                                 paste("median = ",
                                       round(De.stats[i,4], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("median.weighted" %in% summary[j] == TRUE,
                                 paste("weighted median = ",
                                       round(De.stats[i,5], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("kdemax" %in% summary[j] == TRUE,
                                 paste("kdemax = ",
                                       round(De.stats[i,6], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("sdrel" %in% summary[j] == TRUE,
                                 paste("rel. sd = ",
                                       round(De.stats[i,8], 2), " %",
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("sdabs" %in% summary[j] == TRUE,
                                 paste("abs. sd = ",
                                       round(De.stats[i,7], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("serel" %in% summary[j] == TRUE,
                                 paste("rel. se = ",
                                       round(De.stats[i,10], 2), " %",
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("seabs" %in% summary[j] == TRUE,
                                 paste("abs. se = ",
                                       round(De.stats[i,9], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("skewness" %in% summary[j] == TRUE,
                                 paste("skewness = ",
                                       round(De.stats[i,13], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("kurtosis" %in% summary[j] == TRUE,
                                 paste("kurtosis = ",
                                       round(De.stats[i,14], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("in.2s" %in% summary[j] == TRUE,
                                 paste("in 2 sigma = ",
                                       round(sum(data[[i]][,7] > -2 &
                                                   data[[i]][,7] < 2) /
                                               nrow(data[[i]]) * 100 , 1),
                                       " %   ",
                                       sep = ""),
                                 ""),
                          ifelse("sdabs.weighted" %in% summary[j] == TRUE,
                                 paste("abs. weighted sd = ",
                                       round(De.stats[i,15], 2), " %",
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("sdrel.weighted" %in% summary[j] == TRUE,
                                 paste("rel. weighted sd = ",
                                       round(De.stats[i,16], 2), " %",
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("seabs.weighted" %in% summary[j] == TRUE,
                                 paste("abs. weighted se = ",
                                       round(De.stats[i,17], 2), " %",
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("serel.weighted" %in% summary[j] == TRUE,
                                 paste("rel. weighted se = ",
                                       round(De.stats[i,18], 2), " %",
                                       " | ",
                                       sep = ""),
                                 "")
        )
      }

      summary.text <- paste(summary.text, collapse = "")

      label.text[[length(label.text) + 1]]  <- paste(
        "  ",
        summary.text,
        sep = "")
    }

    ## remove outer vertical lines from string
    for(i in 2:length(label.text)) {
      label.text[[i]] <- substr(x = label.text[[i]],
                                start = 3,
                                stop = nchar(label.text[[i]]) - 3)
    }
  }

  ## remove dummy list element
  label.text[[1]] <- NULL

  ## read out additional parameters -------------------------------------------
  if("main" %in% names(list(...))) {
    main <- list(...)$main
  } else {
    main <- expression(bold(paste(D[e], " distribution")))
  }

  if("xlab" %in% names(list(...))) {
    xlab <- list(...)$xlab
  } else {
    xlab <- expression(paste(D[e], " (Gy)"))
  }

  if("ylab" %in% names(list(...))) {
    ylab <- list(...)$ylab
  } else {
    ylab <- c("Density", "Cumulative frequency")
  }

  if("xlim" %in% names(list(...))) {
    xlim.plot <- list(...)$xlim
  } else {
    xlim.plot <- c(min(c(De.global - De.error.global),
                       De.density.range[1],
                       na.rm = TRUE),
                   max(c(De.global + De.error.global),
                       De.density.range[2],
                       na.rm = TRUE))
  }

  if("ylim" %in% names(list(...))) {
    ylim.plot <- list(...)$ylim
  } else {
    ylim.plot <- c(De.density.range[3],
                   De.density.range[4],
                   1,
                   max(De.stats[,1]))
  }

  if("log" %in% names(list(...))) {
    log.option <- list(...)$log
  } else {
    log.option <- ""
  }

  if(length(data) > 1) {
    if("col" %in% names(list(...))) {
      colours <- matrix(rep(list(...)$col, each = 4),
                        nrow = length(data),
                        byrow = TRUE)
    } else {
      colours <- matrix(rep(1:length(data), 4),
                        nrow = length(data))
    }
  } else {
    if("col" %in% names(list(...))) {
      colours <- matrix(c(list(...)$col),
                        nrow = 1)
    } else {
      colours <- matrix(c("#3F489D",
                          "black",
                          "black",
                          "grey90"),
                        nrow = 1)
    }
  }

  if("lty" %in% names(list(...))) {
    lty <- list(...)$lty
  } else {
    lty <- seq(2, 7 * length(data))
  }

  if("cex" %in% names(list(...))) {
    cex <- list(...)$cex
  } else {
    cex <- 1
  }

  if("fun" %in% names(list(...))) {
    fun <- list(...)$fun
  } else {
    fun <- FALSE
  }

  ## convert keywords into summary placement coordinates
  if(missing(summary.pos) == TRUE) {
    summary.pos <- c(xlim.plot[1], ylim.plot[2])
    summary.adj <- c(0, 1)
  } else if(length(summary.pos) == 2) {
    summary.pos <- summary.pos
    summary.adj <- c(0, 1)
  } else if(summary.pos[1] == "topleft") {
    summary.pos <- c(xlim.plot[1], ylim.plot[2])
    summary.adj <- c(0, 1)
  } else if(summary.pos[1] == "top") {
    summary.pos <- c(mean(xlim.plot), ylim.plot[2])
    summary.adj <- c(0.5, 1)
  } else if(summary.pos[1] == "topright") {
    summary.pos <- c(xlim.plot[2], ylim.plot[2])
    summary.adj <- c(1, 1)
  }  else if(summary.pos[1] == "left") {
    summary.pos <- c(xlim.plot[1], mean(ylim.plot[1:2]))
    summary.adj <- c(0, 0.5)
  } else if(summary.pos[1] == "center") {
    summary.pos <- c(mean(xlim.plot), mean(ylim.plot[1:2]))
    summary.adj <- c(0.5, 0.5)
  } else if(summary.pos[1] == "right") {
    summary.pos <- c(xlim.plot[2], mean(ylim.plot[1:2]))
    summary.adj <- c(1, 0.5)
  }else if(summary.pos[1] == "bottomleft") {
    summary.pos <- c(xlim.plot[1], ylim.plot[1])
    summary.adj <- c(0, 0)
  } else if(summary.pos[1] == "bottom") {
    summary.pos <- c(mean(xlim.plot), ylim.plot[1])
    summary.adj <- c(0.5, 0)
  } else if(summary.pos[1] == "bottomright") {
    summary.pos <- c(xlim.plot[2], ylim.plot[1])
    summary.adj <- c(1, 0)
  }




  ## assign polygon coordinates
  polygons <- matrix(nrow = length(data), ncol = 8)

  for(i in 1:length(data)) {
    lims.x <- if(dispersion == "sd") {
      c(De.stats[i,3] - De.stats[i,7],
        De.stats[i,3] - De.stats[i,7],
        De.stats[i,3] + De.stats[i,7],
        De.stats[i,3] + De.stats[i,7])
    } else if(dispersion == "2sd") {
      c(De.stats[i,3] - 2 * De.stats[i,7],
        De.stats[i,3] - 2 * De.stats[i,7],
        De.stats[i,3] + 2 * De.stats[i,7],
        De.stats[i,3] + 2 * De.stats[i,7])
    } else if(dispersion == "qr") {
      c(De.stats[i,11],
        De.stats[i,11],
        De.stats[i,12],
        De.stats[i,12])
    } else {
      rep(NA, 4)
    }

    polygons[i,] <- c(lims.x, c(-2 * ylim.plot[2],
                                2 * ylim.plot[2],
                                2 * ylim.plot[2],
                                -2 * ylim.plot[2]))
  }

  ## plot data sets -----------------------------------------------------------

  ## setup plot area
  if(length(summary) >= 1 & summary.pos[1] == "sub") {
    toplines <- length(data)
  } else {toplines <- 1}

  par(mar = c(4.5, 5.5, 2.5 + toplines, 4.5),
      xpd = FALSE,
      cex = cex)

  ## create empty plot to set plot dimensions
  plot(NA,
       xlim = xlim.plot,
       ylim = ylim.plot[1:2],
       main = "",
       xlab = "",
       ylab = "",
       log = log.option,
       axes = FALSE,
       frame.plot = FALSE)

  ## plot dispersion polygons
  if(length(dispersion) == 1) {
    for(i in 1:length(data)) {
      polygon(x = polygons[i,1:4],
              y = polygons[i,5:8],
              col = polygon.col[i],
              border = FALSE)
    }
  }

  ## plot measures of centrality
  if(length(centrality) >= 1) {
    for(i in 1:length(data)) {
      for(j in 1:length(centrality)) {
        if(centrality[j] == "mean") {
          abline(v = De.stats[i,2], col = colours[i,1], lty = lty[j + 1])
          text(De.stats[i,2] - par()$cxy[1] * 0.5,
               ylim.plot[2] * 0.99, "mean",
               srt = 90, adj = 1, col = colours[i, 2], cex = 0.8 * cex)
        } else if(centrality[j] == "mean.weighted") {
          abline(v = De.stats[i,3], col = colours[i,1], lty = lty[j + 1])
          text(De.stats[i,3] - par()$cxy[1] * 0.5,
               ylim.plot[2] * 0.99, "weighted mean",
               srt = 90, adj = 1, col = colours[i, 2], cex = 0.8 * cex)
        } else if(centrality[j] == "median") {
          abline(v = De.stats[i,4], col = colours[i,1], lty = lty[j + 1])
          text(De.stats[i,4] - par()$cxy[1] * 0.5,
               ylim.plot[2] * 0.99, "median",
               srt = 90, adj = 1, col = colours[i, 2], cex = 0.8 * cex)
        } else if(centrality[j] == "median.weighted") {
          abline(v = De.stats[i,5], col = colours[i,1], lty = lty[j + 1])
          text(De.stats[i,5] - par()$cxy[1] * 0.5,
               ylim.plot[2] * 0.99, "weighted median",
               srt = 90, adj = 1, col = colours[i, 2], cex = 0.8 * cex)
        } else if(centrality[j] == "kdemax") {
          abline(v = De.stats[i,6], col = colours[i,1], lty = lty[j + 1])
          text(De.stats[i,6] - par()$cxy[1] * 0.5,
               ylim.plot[2] * 0.99, "KDE max",
               srt = 90, adj = 1, col = colours[i, 2], cex = 0.8 * cex)
        }
      }
      j <- 1
    }
  }

  ## add probability density plot
  par(new = TRUE)
  plot(NA,
       main     = "",
       xlab     = xlab,
       ylab     = ylab[1],
       xlim     = xlim.plot,
       ylim     = ylim.plot[1:2],
       log      = log.option,
       cex      = cex,
       cex.lab  = cex,
       cex.main = cex,
       cex.axis = cex)

  for(i in 1:length(data)) {
    lines(x = De.density[[i]]$x,
          y = De.density[[i]]$y,
          col = colours[i, 1])
  }

  ## add plot title
  title(main, line = toplines + 1.2)

  ## add summary content
  for(i in 1:length(data)) {
    if(summary.pos[1] != "sub") {
      text(x = summary.pos[1],
           y = summary.pos[2],
           adj = summary.adj,
           labels = label.text[[i]],
           col = colours[i, 2],
           cex = cex * 0.8)
    } else {
      if(mtext == "") {
        mtext(side = 3,
              line = toplines + 0.3 - i,
              text = label.text[[i]],
              col = colours[i, 2],
              cex = cex * 0.8)
      }
    }
  }

  if(values.cumulative == TRUE) {
    ## create empty overlay plot
    par(new = TRUE) # adjust plot options
    plot(NA, # add empty plot, scaled to secondary plot content
         xlim = xlim.plot,
         ylim = ylim.plot[3:4],
         log  = log.option,
         main = "",
         xlab = "",
         ylab = "",
         axes = FALSE,
         frame.plot = FALSE)

    ## add secondary y-axis
    axis(side = 4, labels = TRUE, cex.axis = cex) # add second y-axis
    mtext(ylab[2], side = 4, line = 3, cex = cex) # add second y-axis label

    ## add De error bars
    for(i in 1:length(data)) {
      arrows(data[[i]][,1] - data[[i]][,2]/2,
             1:length(data[[i]][,1]),
             data[[i]][,1] + data[[i]][,2]/2,
             1:length(data[[i]][,1]),
             code = 3,
             angle = 90,
             length = 0.05,
             col = colours[i, 3])

      ## add De measurements
      points(data[[i]][,1], 1:De.stats[i,1],
             col = colours[i, 3],
             pch = 20)
    }
  }

  ## add empty plot
  par(new = TRUE)
  plot(NA,
       ann = FALSE,
       axes = FALSE,
       xlim     = xlim.plot,
       ylim     = ylim.plot[1:2],
       log      = log.option,
       cex      = cex,
       cex.lab  = cex,
       cex.main = cex,
       cex.axis = cex)

  ## FUN by R Luminescence Team
  if(fun==TRUE){sTeve()}

  if(output == TRUE) {
    return(list(De.stats = De.stats,
                summary.pos = summary.pos,
                De.density = De.density))
  }

}
