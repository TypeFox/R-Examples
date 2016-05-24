#' Plot Frequency Profiles.
#' 
#' Use \code{plot_freqprof} to plot frequency profile data generated from 
#' \code{\link{freqprof}}.
#' 
#' @param data.freqprof data formated into class \code{freqprof}.
#' @param yAxis a string labelling the y-axis, defaults to 
#'   \code{data.freqprof$type}.
#' @param xAxisUnits a string indicating x-axis units, defaults to "sec".
#' @param panel.in if \code{FALSE} the first panel of the frequency profile, the
#'   window moving in, is not plotted.
#' @param panel.out if \code{FALSE} the third panel of the frequency profile, 
#'   the window moving out, is not plotted.
#' @param gg if \code{TRUE}, will use \code{ggplot2} to plot frequency profiles.
#' @param multiPlot if \code{TRUE}, will plot each behavior in its own panel.
#' @param tick.every the spacing between each plot tick mark. By default, N/30
#'   where N is the number of time units.
#' @param label.every label every X ticks, where X = label.every. By default, 
#'   label.every = 3.
#' @return Returns a frequency profiles plot.
#' @importFrom reshape2 melt
#' @importFrom graphics abline axis lines mtext par plot
#' @import ggplot2
#' @export
#' @examples
#' data(s58)
#' plot_freqprof(freqprof(s58))
plot_freqprof = function(data.freqprof,
                         yAxis      = NULL,
                         xAxisUnits = "sec",
                         panel.in   = TRUE,
                         panel.out  = TRUE,
                         gg         = FALSE,
                         multiPlot  = FALSE,
                         tick.every = round(length(data.freqprof$data$time) / 31),
                         label.every = 3) {
  
  # Extract relevant data from data.freqprof
  res      <- data.freqprof$data
  panels   <- res$panels
  res      <- res[, -2]
  freqprof <- res[, -1]
  t        <- res$time
  
  window     <- data.freqprof$window
  step       <- data.freqprof$step
  resolution <- data.freqprof$resolution
  type       <- data.freqprof$type
  
  # Panels limits
  x.panel.left  <- max(which(data.freqprof$data$panels == 1)) * resolution
  x.panel.right <- min(which(data.freqprof$data$panels == 3)) * resolution
  
  # X-axis limits
  xmin <- ifelse(test = panel.in,  yes = min(t), no = x.panel.left)
  xmax <- ifelse(test = panel.out, yes = max(t), no = x.panel.right)
  
  # If no custom yAxis label given, label according to data.freqprof$type
  if(is.null(yAxis)) {
    yAxis = switch(type,
                   sum        = 'Moving sum',
                   proportion = 'Moving proportion')
  }
  
  # Color-blind friendly palette
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  if(gg) {

    res.melt <- melt(res, id = "time")
    
    # Graphing function
    p <- ggplot_fp(data1       = res.melt, 
                   resolution  = resolution, 
                   step        = step,
                   yAxis       = yAxis,
                   xAxisUnits  = xAxisUnits,
                   xmin        = xmin,
                   xmax        = xmax,
                   tick.every  = tick.every,
                   label.every = label.every)
    
    if(panel.in) {
      p = p + geom_vline(xintercept = x.panel.left)
    }
    
    if(panel.out) {
      p = p + geom_vline(xintercept = x.panel.right) 
    }
    
    if (multiPlot) {
      p = p + facet_grid(variable ~ .) + theme(legend.position = "none")
    }
    
    print(p)
  } else {
    
    # no ggplot
    
    if(is.null(ncol(freqprof))) {
      # case of only one column selected
      freqprof = as.data.frame(freqprof)
    }
    
    if(multiPlot) {
      plotBehavior = function(j) {
        plot(x    = t,
             y    = freqprof[, j],
             type = 'l',
             col  = cbbPalette[(j %% length(cbbPalette)) + 1],
             ylim = c(0, max(freqprof)),
             xlim = c(xmin, xmax),
             xlab = '',
             ylab = '',
             xaxt = 'n')
        n.minor = length(ax <- seq(from = xmin, 
                                   to   = xmax, 
                                   by   = tick.every))
        axis(1, at = ax, labels = FALSE)
        lab = ax[seq(1, n.minor, by = label.every)]
        axis(1, at = lab, labels = TRUE)
        if(panel.in) {
          abline(v = x.panel.left)
        } 
        if(panel.out) {
          abline(v = x.panel.right)
        }
      }
      
      par(mfrow = c(ncol(freqprof), 1))
      par(cex   = 0.6)
      par(mar   = c(2, 2, 2, 2), 
          oma   = c(4, 4, 0.5, 0.5))
      
      for(j in 1:ncol(freqprof)) plotBehavior(j)
      
      mtext(paste('Time (', resolution * step, ' ', xAxisUnits,')', 
                  sep = ""), 
            side  = 1, 
            outer = TRUE, 
            cex   = 0.7, 
            line  = 2.2)
      mtext(yAxis, 
            side  = 2, 
            outer = TRUE, 
            cex   = 0.7, 
            line  = 2.2)
    }
    else {
      plot(x    = t,
           y    = freqprof[, 1],
           type = 'l',
           col  = cbbPalette[1],
           ylim = c(0, max(freqprof)),
           xlim = c(xmin, xmax),
           ylab = yAxis,
           xlab = paste('Time (', resolution, ' ', xAxisUnits, ')', sep = ""),
           xaxt = 'n')
      n.minor = length(ax <- seq(from = xmin, 
                                 to   = xmax, 
                                 by   = tick.every))
      axis(1, at = ax, labels = FALSE)
      lab = ax[seq(1, n.minor, by = label.every)]
      axis(1, at = lab, labels = TRUE)
      
      if(ncol(freqprof) > 1) {
        for(j in 2:ncol(freqprof)) {
          lines(t,
                freqprof[, j],
                type = 'l',
                col = cbbPalette[(j %% length(cbbPalette)) + 1])
        }
      }
      if(panel.in) {
        abline(v = x.panel.left)
      } 
      if(panel.out) {
        abline(v = x.panel.right)
      }
    }
    
  }
}

#' Internal ggplot Wrapper to Graph Frequency Profiles
#' 
#' @param data1 data formated into \code{freqprof} class.
#' @param resolution resolution of \code{freqprof} data
#' @param step step size of \code{freqprof} data
#' @param yAxis a string providing a label for the y-axis.
#' @param xAxisUnits a string indicating which unit has been used. By default, 
#'   "sec".
#' @param xmin x-axis minimum value
#' @param xmax x-axis maximum value
#' @param tick.every the spacing between each tick. By default, N/30 where N is 
#'   the number of time units.
#' @param label.every label every X ticks, where X = label.every. By default, 
#'   label.every = 3.
#' @return A ggplot of the frequency profile data in \code{data1}
#' 
ggplot_fp <- function(data1, 
                      resolution  = resolution, 
                      step        = step,
                      yAxis       = yAxis,
                      xAxisUnits  = xAxisUnits,
                      xmin        = xmin,
                      xmax        = xmax,
                      tick.every  = tick.every,
                      label.every = label.every) {
  
  p <- with(data1, {
     ggplot(data1,
            aes(x      = time, 
                y      = value,
                color  = variable, 
                group  = variable)) +
      geom_line(size = 0.8) +
      labs(title = "Frequency Profile") +
      xlab(paste('Time (', resolution * step, ' ', xAxisUnits, ')', sep = "")) + 
      ylab(paste(yAxis)) +
      scale_x_continuous(limits       = c(xmin, xmax),
                        minor_breaks = round(seq(xmin, xmax, by = tick.every)),
                        breaks       = round(seq(xmin, xmax,
                                             by = tick.every * label.every))) +
      scale_color_discrete(name = "Behavior") +
      theme(title = element_text(size = 17, face = "bold"),
            axis.text.x  = element_text(size   = 12,
                                        color  = "#3f3f3f",
                                        margin = margin(t = 0.4, unit = "cm")),
            axis.text.y  = element_text(size   = 12, 
                                        color  = "#3f3f3f",
                                        margin = margin(r = 0.4, unit = "cm")),
            axis.title.x = element_text(size   = 14, 
                                        face   = "bold",
                                        margin = margin(t = 0.4, unit = "cm")),
            axis.title.y = element_text(size   = 14, 
                                        face   = "bold",
                                        margin = margin(r = 0.4, unit = "cm")),
            legend.text       = element_text(size   = 12),
            panel.background  = element_rect(fill   = "#f6f6f6"),
            panel.grid.major  = element_line(color  = "#e9e9e9"),
            panel.grid.minor  = element_line(color  = "#e9e9e9"),
            axis.line         = element_line(color  = "#a8a8a8"),
            axis.ticks        = element_line(color  = "black", size = 0.5),
            axis.ticks.length = unit(-0.2, "cm"))
  })
  return(p)
}