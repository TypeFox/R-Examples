#' Plotting method for \code{reproData} objects
#'
#' Plots the cumulated number of offspring as a
#' function of either time and concentration, time only (for a fixed
#' concentration), concentration only (for a given target time). If both
#' concentration and target time are fixed, the function additionally plots
#' the experimental values for the minimum available concentration.
#'
#' @param x an object of class \code{reproData}
#' @param xlab a title for the \eqn{x}-axis (optional)
#' @param ylab a title for the \eqn{y}-axis
#' @param main main title for the plot
#' @param target.time a numeric value corresponding to some observed time in \code{data}
#' @param concentration a numeric value corresponding to some concentration in \code{data}
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param pool.replicate if \code{TRUE}, the datapoints of each replicate are
#' summed for a same concentration
#' @param log.scale if \code{TRUE}, displays \eqn{x}-axis in log-scale
#' @param addlegend if \code{TRUE}, adds a default legend to the plot
#' @param remove.someLabels if \code{TRUE}, removes 3/4 of X-axis labels in
#' \code{'ggplot'} style to avoid the label overlap
#' @param \dots Further arguments to be passed to generic methods.
#' @note When \code{style = "ggplot"}, the function calls package
#' \code{\link[ggplot2]{ggplot}} and returns an object of class \code{ggplot}.
#' @keywords plot
#'
#' @examples
#'
#' library(ggplot2)
#'
#' # (1) Load the data
#' data(cadmium1)
#' cadmium1 <- reproData(cadmium1)
#'
#' # (2) Plot the reproduction data
#' plot(cadmium1)
#'
#' # (3) Plot the reproduction data for a fixed time with a ggplot style
#' plot(cadmium1, target.time = 21, style = "ggplot")
#'
#' # (4) Plot the reproduction data for a fixed concentration
#' plot(cadmium1, concentration = 4.36, style = "ggplot")
#'
#' # (5) Plot the reproduction data for a fixed concentration and target.time
#' plot(cadmium1, target.time = 21, concentration = 0.86)
#'
#' @import ggplot2
#' @import grDevices
#' @importFrom graphics plot axis legend lines par points polygon
#' segments title
#' @importFrom methods is
#' @importFrom stats aggregate
#' 
#' @export
plot.reproData <- function(x,
                           xlab,
                           ylab = "Cumulated Number of offspring",
                           main = NULL,
                           target.time = NULL,
                           concentration = NULL,
                           style = "generic",
                           pool.replicate = FALSE,
                           log.scale = FALSE,
                           addlegend = FALSE,
                           remove.someLabels = FALSE, ...) {
  if(! is(x, "reproData"))
    stop("plot.reproData: object of class reproData expected")

  if (pool.replicate) {
    # agregate by sum of replicate
    x <- cbind(aggregate(cbind(Nreprocumul, Nsurv, Ninit) ~ time + conc, x, sum),
               replicate = 1)
  }

  if (is.null(target.time) && is.null(concentration)) {
    reproDataPlotFull(x, xlab, ylab, style, remove.someLabels)
  }
  else if (! is.null(target.time) && is.null(concentration)) {
    reproDataPlotTargetTime(x, xlab, ylab, main, target.time,
                            style, log.scale, addlegend,
                            remove.someLabels)
  }
  else if (is.null(target.time) && ! is.null(concentration)) {
    reproDataPlotFixedConc(x, xlab, ylab, main, concentration, style, addlegend,
                           remove.someLabels)
  }
  else {
    reproDataPlotReplicates(x, xlab, ylab, target.time, concentration,
                            style, addlegend)
  }
}



reproDataPlotFull <- function(data, xlab, ylab, style = "generic",
                              remove.someLabels) {
  dataPlotFull(data, xlab, ylab, "Nreprocumul", style,
               remove.someLabels)
}


#' @import ggplot2
#' @importFrom dplyr filter
reproDataPlotTargetTime <- function(x,
                                    xlab,
                                    ylab,
                                    main,
                                    target.time,
                                    style,
                                    log.scale,
                                    addlegend,
                                    remove.someLabels) {

  if (missing(xlab)) xlab <-"Concentration"

  # plot of cumulated number of offspring as a funtion of concentration
  # for a fixed time

  if (!target.time %in% x$time)
    stop("[target.time] is not one of the possible time !")

  # select the target.time
  xf <- filter(x, x$time == target.time)

  # Selection of datapoints that can be displayed given the type of scale
  sel <- if(log.scale) xf$conc > 0 else TRUE
  x <- xf[sel, ]
  transf_data_conc <- optLogTransform(log.scale, x$conc)

  # Concentration values used for display in linear scale
  display.conc <- (function() {
    x <- optLogTransform(log.scale, x$conc)
    if(log.scale) exp(x) else x
  })()

  # Define visual parameters
  mortality <- c(0, 1) # code 0/1 mortality
  nomortality <- match(x$Nsurv == x$Ninit, c(TRUE, FALSE))

  # without mortality
  mortality <- mortality[nomortality] # vector of 0 and 1

  # encodes mortality empty dots (1) and not mortality solid dots (19)
  if (style == "generic") {
    mortality[which(mortality == 0)] <- 19
  }
  if (style == "ggplot") {
    mortality[which(mortality == 0)] <- "No"
    mortality[which(mortality == 1)] <- "Yes"
  }

  # default legend argument
  legend.position <- "right"
  legend.title <- "Mortality"
  legend.name.no <- "No"
  legend.name.yes <- "Yes"

  # generic
  if (style == "generic") {
    plot(transf_data_conc,
         x$Nreprocumul,
         xlab = xlab,
         ylab = ylab,
         main = main,
         pch = mortality,
         yaxt = "n",
         xaxt = "n")
    # axis
    axis(side = 2, at = pretty(c(0, max(x$Nreprocumul))))
    axis(side = 1, at = transf_data_conc,
         labels = display.conc)

    # legend
    if (addlegend) {
      legend(legend.position,title = legend.title, pch = c(19, 1), bty = "n",
             legend = c(legend.name.no, legend.name.yes))
    }
  }

  #ggplot2
  if (style == "ggplot") {
    df <- data.frame(x,
                     transf_data_conc,
                     display.conc,
                     Mortality = mortality)

    # plot
    gp <- ggplot(df, aes(transf_data_conc, Nreprocumul,
                         fill = Mortality)) +
      geom_point(size = 3, pch = 21) +
      scale_fill_manual(values = if (
        length(levels(df$Mortality)) == 1 && levels(df$Mortality) == "Yes") {
        "white"
        } else if (length(levels(df$Mortality)) == 1 &&
                   levels(df$Mortality) == "No") {
          "black"
        } else {
          c("black", "white")
        }) +
      labs(x = xlab, y = ylab) +
      ggtitle(main) +
      theme_minimal() +
      scale_colour_hue(legend.title, breaks = c("No","Yes"),
                       labels = c(legend.name.no, legend.name.yes)) +
      scale_x_continuous(breaks = df$transf_data_conc,
                         labels = if (remove.someLabels) {
                           exclude_labels(df$display.conc)
                         } else {
                           df$display.conc
                         })

    if (addlegend) {
      return(gp)
    } else {
      gp <- gp + theme(legend.position = "none")
      return(gp)
    }
  }
}

reproDataPlotFixedConc <- function(x,
                                   xlab,
                                   ylab,
                                   main,
                                   concentration,
                                   style = "generic",
                                   addlegend = FALSE,
                                   remove.someLabels = FALSE) {
  dataPlotFixedConc(x, xlab, ylab, main, "Nreprocumul",
                    concentration, style, addlegend, remove.someLabels)
}

reproDataPlotReplicates <- function(x,
                                    xlab,
                                    ylab,
                                    target.time,
                                    concentration,
                                    style,
                                    addlegend) {
  dataPlotReplicates(x, xlab, ylab, "Nreprocumul", target.time,
                     concentration, style, addlegend)
}
