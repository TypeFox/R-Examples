#' Plotting method for \code{survData} objects
#'
#' Plots the number of survivors as a
#' function of either time and concentration, time only (for a fixed
#' concentration), concentration only (for a given target time). If both
#' concentration and target time are fixed, the function additionally plots
#' the experimental values for the minimum available concentration.
#'
#' @param x an object of class \code{survData}
#' @param xlab a title for the \eqn{x}-axis (optional)
#' @param ylab a label for the \eqn{y}-axis
#' @param main main title for the plot
#' @param target.time a numeric value corresponding to some observed time in \code{data}
#' @param concentration a numeric value corresponding to some concentration in \code{data}
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param pool.replicate if \code{TRUE}, the datapoints of each replicate are
#' summed for a same concentration
#' @param log.scale if \code{TRUE}, displays \eqn{x}-axis in log scale
#' @param addlegend if \code{TRUE}, adds a default legend to the plot
#' @param remove.someLabels if \code{TRUE}, removes 3/4 of X-axis labels in
#' \code{'ggplot'} style to avoid label overlap
#' @param \dots Further arguments to be passed to generic methods.
#' @note When \code{style = "ggplot"}, the function calls package
#' \code{\link[ggplot2]{ggplot}} and returns an object of class \code{ggplot}.
#'
#' @keywords plot
#'
#' @examples
#'
#' library(ggplot2)
#'
#' # (1) Load the data
#' data(zinc)
#' zinc <- survData(zinc)
#'
#' # (2) Plot survival data
#' plot(zinc)
#'
#' # (3) Plot survival data with a ggplot style
#' plot(zinc, style = "ggplot", remove.someLabels = TRUE)
#' 
#' # (4) Plot the survival data for a fixed time with a ggplot style
#' plot(zinc, target.time = 21, style = "ggplot")
#'
#' # (5) Plot the survival data for a fixed concentration
#' plot(zinc, concentration = 0.66, addlegend = TRUE)
#'
#' # (6) Plot survival data for a fixed concentration and
#' # target.time with ggplot style
#' plot(zinc, style = "ggplot", target.time = 21, concentration = 0.66)
#'
#' @import ggplot2
#' @import grDevices
#' @importFrom graphics plot axis legend lines par points polygon
#' segments title
#' @importFrom methods is
#' @importFrom stats aggregate
#'
#' @export
plot.survData <- function(x,
                          xlab,
                          ylab = "Number of surviving individuals",
                          main = NULL,
                          target.time = NULL,
                          concentration = NULL,
                          style = "generic",
                          pool.replicate = FALSE,
                          log.scale = FALSE,
                          addlegend = FALSE,
                          remove.someLabels = FALSE, ...) {

  if(! is(x,"survData"))
    stop("plot.survData: object of class survData expected")

  if (pool.replicate) {
    # agregate by sum of replicate
    x <- cbind(aggregate(Nsurv ~ time + conc, x, sum),
               replicate = 1)
  }

  if (is.null(target.time) && is.null(concentration)) {
    survDataPlotFull(x, xlab, ylab, style, remove.someLabels)
  }
  else if (! is.null(target.time) && is.null(concentration)) {
    survDataPlotTargetTime(x, xlab, ylab, main, target.time,
                           style, log.scale, addlegend, remove.someLabels)
  }
  else if (is.null(target.time) && ! is.null(concentration)) {
    survDataPlotFixedConc(x, xlab, ylab, main, concentration,
                          style, addlegend, remove.someLabels)
  }
  else {
    survDataPlotReplicates(x, xlab, ylab, target.time, concentration, style,
                           addlegend)
  }
}


# [ReplicateIndex(data)] builds a list of indices, each one named after
# a replicate of [data], thus providing a dictionary from replicate names to
# integer keys.
ReplicateIndex <- function(data) {
  replicate <- unique(data$replicate)
  r <- as.list(seq(1, length(replicate)))
  names(r) <- as.character(replicate)
  return(r)
}


# [plotMatrixGeometry(n)] returns a vector [c(w,h)] such that a matrix of plots
# of dimension ([w], [h]) is big enough to display [n] plots in a pretty way.
# This will typically be used in [par(mfrow)] calls.
plotMatrixGeometry <- function(nblevels) {
  PlotPar <- c(c(2, 2), c(2, 3), c(2, 4), c(3, 3), c(2, 5), c(3, 4), c(3, 5),
               c(4, 4))
  NbPlotTheo <- matrix(ncol = 2, nrow = 8)
  NbPlotTheo[, 1] <- c(1, 3, 5, 7, 9, 11, 13, 15)
  NbPlotTheo[, 2] <- c(4, 6, 8, 9, 10, 12, 15, 16)
  if (nblevels < 15) {
    i <- NbPlotTheo[NbPlotTheo[, 2] - nblevels > 0, 1][1]
  } else {
    i <- 15
  }
  return(c(PlotPar[i], PlotPar[i + 1]))
}

# General full plot: one subplot for each concentration, and one color for
# each replicate (for generic graphics)
dataPlotFullGeneric <- function(data, xlab, ylab, resp) {
  replicate.index <- ReplicateIndex(data)

  # creation of a vector of colors
  colors <- rainbow(length(unique(data$replicate)))
  pchs <- as.numeric(unique(data$replicate))
  # split of the graphical window in subplots
  par(mfrow = plotMatrixGeometry(length(unique(data$conc))))

  by(data, data$conc, function(x) {
    # bakground
    plot(x$time, rep(0, length(x$time)),
         xlab = xlab,
         ylab = ylab,
         ylim = c(0, max(x[, resp])),
         type = "n",
         col = 'white',
         xaxt = "n",
         yaxt = "n")

    # axis
    axis(side = 1, at = sort(unique(x[, "time"])))
    axis(side = 2, at = unique(round(pretty(c(0, max(x[, resp]))))))

    # lines and points
    by(x, x$replicate, function(y) {
      index <- replicate.index[[y$replicate[1]]]
      lines(y$time, y[, resp],
            type = "l",
            col = colors[index])
      points(y$time, y[, resp],
             pch = pchs[index],
             col = colors[index])
    })

    # title
    title(paste("Conc: ", unique(x$conc), sep = ""))
  })
  
  par(mfrow = c(1, 1))
}

# general full plot (ggplot variant): one subplot for each concentration,
# and one color for each replicate
#' @import ggplot2
dataPlotFullGG <- function(data, xlab, ylab, resp, remove.someLabels) {
  
  time = NULL
  Nsurv = NULL
  
  data$response <- data[,resp]
  
  # create ggplot object Nsurv / time / replicate / conc
  fg <- ggplot(data, aes(time, response, colour = factor(replicate))) +
    geom_point() +
    geom_line() +
    labs(x = xlab, y = ylab) +
    facet_wrap(~conc, ncol = 2) +
    scale_x_continuous(breaks = unique(data$time),
                       labels = if (remove.someLabels) {
                         exclude_labels(unique(data$time))
                       } else {
                         unique(data$time)
                       }
    ) +
    scale_y_continuous(breaks = unique(round(pretty(c(0, max(data[, resp])))))) +
    expand_limits(x = 0, y = 0) +
    theme_minimal()
  
  fd <- fg + theme(legend.position = "none") # remove legend
  
  return(fd)
}

dataPlotFull <- function(data, xlab, ylab, resp, style = "generic",
                         remove.someLabels = FALSE) {

  if (missing(xlab)) xlab <- "Time"

  if (style == "generic")
    dataPlotFullGeneric(data, xlab, ylab, resp)
  else if (style == "ggplot")
    dataPlotFullGG(data, xlab, ylab, resp, remove.someLabels)
  else stop("Unknown plot style")
}

survDataPlotFull <- function(data, xlab, ylab, style = "generic",
                             remove.someLabels = FALSE) {
  dataPlotFull(data, xlab, ylab, "Nsurv", style, remove.someLabels)
}

#' @import ggplot2
#' @importFrom dplyr %>% filter
survDataPlotTargetTime <- function(x, xlab, ylab, main, target.time,
                                   style, log.scale, addlegend,
                                   remove.someLabels) {
  if (missing(xlab)) xlab <-"Concentration"

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

  # vector color
  x$color <- as.numeric(as.factor(x$replicate))

  if (style == "generic") {
    plot(transf_data_conc, seq(0, max(x$Nsurv),
                               length.out = length(transf_data_conc)),
         type = "n",
         xaxt = "n",
         yaxt = "n",
         main = main,
         xlab = xlab,
         ylab = ylab)

    axis(side = 1, at = transf_data_conc,
         labels = display.conc)
    axis(side = 2, at = unique(round(pretty(c(0, max(x$Nsurv))))),
         labels = unique(round(pretty(c(0, max(x$Nsurv))))))

    # points
    if (length(unique(x$replicate)) == 1) {
      # points
      points(transf_data_conc, x$Nsurv,
             pch = 16)
    } else {
      tt <- xyTable(transf_data_conc, x$Nsurv)
      points(tt$x, tt$y,
             cex = (tt$number) / 3,
             pch = 16)
      if (addlegend) {
        legend("bottomleft",
               legend = sort(unique(tt$number)),
               pt.cex = sort(unique((tt$number) / 3)),
               title = "Number of replicates",
               pch = 16,
               bty = "n")
      }
    }
  }
  else if (style == "ggplot") {
    df <- data.frame(x,
                     transf_data_conc,
                     display.conc)

    if (length(unique(df$replicate)) == 1) {
      gp <- ggplot(df, aes(x = transf_data_conc, y = Nsurv))
    } else {
      gp <- ggplot(df, aes(x = transf_data_conc, y = Nsurv)) +
        stat_sum(aes(size = factor(..n..))) +
        scale_size_discrete("Number of replicates")
    }
    fd <- gp + geom_point() + ggtitle(main) +
      theme_minimal() +
      labs(x = xlab,
           y = ylab) +
      scale_x_continuous(breaks = df$transf_data_conc,
                         labels = if (remove.someLabels) {
                           exclude_labels(df$display.conc)
                         } else {
                           df$display.conc
                         }
      ) +
      scale_y_continuous(breaks = unique(round(pretty(c(0, max(df$Nsurv)))))) +
      expand_limits(x = 0, y = 0)

    # legend option
    if (addlegend) {
      fd
    } else {
      fd + theme(legend.position = "none") # remove legend
    }
  }
  else stop("Unknown plot style")
}

dataPlotFixedConc <- function(x,
                              xlab,
                              ylab,
                              main,
                              resp,
                              concentration,
                              style = "generic",
                              addlegend = FALSE,
                              remove.someLabels = FALSE) {

  if (missing(xlab)) xlab <- "Time"

  legend.position <- ifelse(resp == "Nsurv", "bottomleft", "topleft")

  # check concentration value
  if (!concentration %in% x$conc)
    stop("The argument [concentration] should correspond to one of the tested concentrations")

  # select the concentration
  x <- filter(x, x$conc == concentration)

  # vector color
  x$color <- as.numeric(as.factor(x$replicate))

  if (style == "generic") {
    plot(x$time, x[, resp],
         type = "n",
         xaxt = "n",
         yaxt = "n",
         main = main,
         xlim = range(x$time),
         ylim = c(0, max(x[, resp])),
         xlab = xlab,
         ylab = ylab)

    # one line by replicate
    by(x, list(x$replicate),
       function(x) {
         lines(x$time, x[,resp], # lines
               col = x$color)
         points(x$time, x[,resp], # points
                pch = 16,
                col = x$color)
       })

    # axis
    axis(side = 1, at = sort(unique(x[, "time"])))
    axis(side = 2, at = unique(round(pretty(c(0, max(x[, resp]))))))

    if (addlegend && !unique(x$replicate) == 1) {
      legend(legend.position, legend = unique(x$replicate) ,
             col = unique(x$color),
             pch = 16,
             lty = 1)
    }
  }
  else if (style == "ggplot") {
    x$response <- x[,resp]

    if (length(unique(x$replicate)) == 1) {
      df <- ggplot(x, aes(x = time, y = response))
    } else {
      df <- ggplot(x, aes(x = time, y = response,
                          color = factor(replicate),
                          group = replicate))
    }
    fd <- df + geom_line() + geom_point() + ggtitle(main) +
      theme_minimal() +
      labs(x = xlab,
           y = ylab) +
      scale_color_hue("Replicate") +
      scale_x_continuous(breaks = unique(x$time),
                         labels = if (remove.someLabels) {
                           exclude_labels(unique(x$time))
                         } else {
                           unique(x$time)
                         }) +
      scale_y_continuous(breaks = unique(round(pretty(c(0, max(x$response)))))) +
      expand_limits(x = 0, y = 0)

    if (addlegend) {# only if pool.replicate == FALSE
      fd
    } else {
      fd + theme(legend.position = "none") # remove legend
    }
  }
  else stop("Unknown plot style")
}

survDataPlotFixedConc <- function(x,
                                  xlab,
                                  ylab,
                                  main,
                                  concentration,
                                  style = "generic",
                                  addlegend = FALSE,
                                  remove.someLabels = FALSE) {

  dataPlotFixedConc(x, xlab, ylab, main, "Nsurv", concentration,
                    style, addlegend, remove.someLabels)
}

#' @importFrom dplyr %>% filter
dataPlotReplicates <- function(x,
                               xlab,
                               ylab,
                               resp,
                               target.time,
                               concentration,
                               style,
                               addlegend) {

  if (missing(xlab)) xlab <- "Replicate"

  # check [target.time] and [concentration]
  if (!target.time %in% x$time)
    stop("The argument [target.time] should correspond to one of the observed time points")

  if (!concentration %in% x$conc)
    stop("The argument [concentration] should correspond to one of the tested concentrations")

  # select for concentration and target.time
  xtt <- filter(x, conc == concentration & time == target.time)
  control <- filter(x, conc == min(x$conc) & time == target.time)

  if (style == "generic") {
    par(mfrow = c(1, ifelse(concentration == 0, 1, 2)))

    plot(as.numeric(control$replicate), control[,resp],
         xlab = xlab,
         ylab = ylab,
         main = "Control",
         pch = 16,
         ylim = c(0, max(control[, resp])),
         xaxt = "n",
         yaxt = "n")

    # axis
    axis(side = 1, at = sort(unique(control[, "replicate"])))
    axis(side = 2, at = unique(round(pretty(c(0, max(control[, resp]))))))

    # fixed concentration
    if (! concentration == 0) {
      plot(as.numeric(xtt$replicate), xtt[,resp],
           xlab = xlab,
           ylab = ylab,
           main = paste("Concentration: ", concentration, sep = ""),
           pch = 16,
           ylim = c(0, max(control[, resp])),
           xaxt = "n",
           yaxt = "n")

      # axis
      axis(side = 1, at = sort(unique(xtt[, "replicate"])))
      axis(side = 2, at = unique(round(pretty(c(0, max(xtt[, resp]))))))
    }
  }

  else if (style == "ggplot") {
    dataall <- rbind(control, xtt)
    dataall$response <- dataall[,resp]
    df <- ggplot(dataall, aes(x = replicate, y = response))
    df + geom_point() + labs(x = xlab, y = ylab) +
      scale_x_discrete(breaks = dataall$replicate,
                       labels = dataall$replicate) +
      scale_y_discrete(breaks = unique(round(pretty(c(0, max(dataall$response)))))) +
      expand_limits(x = 0, y = 0) +
      facet_wrap(~conc) + theme_minimal()
  }
  else stop("Unknown plot style")
}

survDataPlotReplicates <- function(x,
                                   xlab,
                                   ylab,
                                   target.time,
                                   concentration,
                                   style,
                                   addlegend) {
  dataPlotReplicates(x, xlab, ylab, "Nsurv", target.time, concentration, style,
                     addlegend)
}
