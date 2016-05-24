#' Plot the curves
#'
#' \code{plotcurves_} is the standard evaluation SE function associated
#' to the non-standard evaluation NSE function \code{plotcurves}.
#' \href{http://adv-r.had.co.nz/Computing-on-the-language.html}{SE functions can be more easily called from other functions.}
#' In SE functions, you need to quote the names of the variables.
#' @param qp output from quickpsy
#' @param panel Name of the variable to be split in panels.
#' @param xpanel Name of the variable to be split in horizontal panels.
#' @param ypanel Name of the variable to be split in vertical panels.
#' @param color Name of the variable codded by color.
#' @param averages If \code{FALSE} averaged probabilities are not plotted
#' (default is \code{TRUE}).
#' @param curves If \code{FALSE} curves are not plotted
#' (default is \code{TRUE})
#' @param thresholds If \code{FALSE} thresholds  are not plotted
#' (default is \code{TRUE})
#' @param ci If \code{FALSE} confidence intervals are not plotted
#' (default is \code{TRUE})
#' @seealso \code{\link{plotcurves}}
#' @examples
#' library(MPDiR) # contains the Vernier data
#' data(Vernier) # ?Venier for the reference
#' fit <- quickpsy(Vernier, Phaseshift, NumUpward, N,
#'                 grouping = .(Direction, WaveForm, TempFreq), B = 10)
#'
#' plotcurves_(fit, xpanel = 'Direction')
#' plotcurves_(fit, color = 'Direction')
#' plotcurves_(fit, xpanel = 'Direction', color = 'WaveForm', ci = FALSE)
#' @import ggplot2
#' @export
plotcurves_ <- function(qp, panel = NULL, xpanel = NULL, ypanel = NULL,
               color = NULL, averages = TRUE, curves = TRUE,
               thresholds = TRUE, ci = TRUE) {

  if (!('thresholds' %in% names(qp))) thresholds <- FALSE
  if (!('thresholdsci' %in% names(qp))) ci <- FALSE

  if (is.logical(qp$guess)) qp$guess <- 0
  if (is.logical(qp$lapses)) qp$lapses <- 0

  p <- ggplot()

  if (qp$log) {
    xmin <- min(qp$averages[[qp$x]])
    xmax <- max(qp$averages[[qp$x]])
    breaks <- signif(exp( seq(log(xmin), log(xmax), length.out=4) ), digits = 2)
    p <- p + scale_x_log10(breaks = breaks)
  }

  groups <- qp$groups
  ngroup <- length(groups)

  if (ngroup == 1) { ###########################################################
    if (!is.null(color)) groups <- setdiff(groups, color)
    if (!is.null(xpanel)) groups <- setdiff(groups, xpanel)
    if (!is.null(ypanel)) groups <- setdiff(groups, ypanel)
    if (!is.null(panel)) groups <- setdiff(groups, panel)
    if (length(groups) == 1) color <- groups[[1]]

    if (!is.null(panel)) p <- p +
      facet_wrap(as.formula(paste0('~',panel)))
    if (!is.null(xpanel)) p <- p +
      facet_grid(as.formula(paste0('.~',xpanel)))
    if (!is.null(ypanel)) p <- p +
      facet_grid(as.formula(paste0(ypanel,'~.')))
  }

  if (ngroup == 2) { ###########################################################
   if (!is.null(color)) groups <- setdiff(groups, color)
   if (!is.null(xpanel)) groups <- setdiff(groups, xpanel)
   if (!is.null(ypanel)) groups <- setdiff(groups, ypanel)
   if (!is.null(panel)) groups <- setdiff(groups, panel)
   if (is.null(color) && length(groups) >= 1) {
     color <- groups[[1]]
     groups <- setdiff(groups,groups[[1]])
   }

   if (is.null(xpanel) && is.null(ypanel) && is.null(panel)) {
     panel <- groups[[1]]
     p <- p + facet_wrap(as.formula(paste0('~',panel)))
   }
   else {
     if (!is.null(panel)) p <- p +
       facet_wrap(as.formula(paste0('~',panel)))
     if (!is.null(xpanel)) p <- p +
       facet_grid(as.formula(paste0('.~',xpanel)))
     if (!is.null(ypanel)) p <- p +
       facet_grid(as.formula(paste0(ypanel,'~.')))
     if (!is.null(xpanel) && !is.null(ypanel)) p <- p +
       facet_grid(as.formula(paste0(ypanel,'~', xpanel)))
   }
  }

  if (ngroup == 3) { ###########################################################
   if (!is.null(color)) groups <- setdiff(groups, color)
   if (!is.null(xpanel)) groups <- setdiff(groups, xpanel)
   if (!is.null(ypanel)) groups <- setdiff(groups, ypanel)
   if (is.null(color)) {
     color <- groups[[1]]
     groups <- setdiff(groups,groups[[1]])
   }
   if (is.null(xpanel)) {
     xpanel <- groups[[1]]
     groups <- setdiff(groups,groups[[1]])
   }
   if (is.null(ypanel)) {
     ypanel <- groups[[1]]
   }

   p <- p + facet_grid(as.formula(paste0(ypanel,'~',xpanel)))
  }

### plotting ###################################################################
  if (ngroup == 0) {
   if (averages) p <- p + geom_point(data = qp$averages,
                          aes_string(x = qp$x, y = 'prob'))
   if (curves) p <- p + geom_line(data = qp$curves,
                        aes_string(x = 'x', y = 'y'))
   if (thresholds) p <- p + geom_linerange(data = qp$thresholds,
                        aes_string(x = 'thre', ymin = qp$guess,
                            ymax = qp$thresholds$prob))
   if (ci) p <- p + geom_errorbarh(data = qp$thresholdsci,
                    height = .03, aes_string(x = 'threinf', xmin = 'threinf',
                    xmax = 'thresup', y = qp$thresholds$prob))
  }
  if (ngroup == 1 || ngroup ==2 || ngroup == 3) {
    if (!is.null(color)) {
      qp$averages[[color]] <- factor(qp$averages[[color]])
      qp$curves[[color]] <- factor(qp$curves[[color]])

      if (averages) p <- p + geom_point(data = qp$averages,
                  aes_string(x = qp$x, y = 'prob', color = color))

      if (curves) p <- p + geom_line(data = qp$curves,
                    aes_string(x = 'x', y = 'y', color = color))

      if (thresholds) {
        qp$thresholds[[color]] <- factor(qp$thresholds[[color]])
        # get present axis limits
        axisYrange <- ggplot_build(p)$panel$ranges[[1]]$y.range
        p <- p + geom_linerange(data = qp$thresholds,
                    aes_string(x = 'thre',
                    		   ymin = axisYrange[1] - .2, #make sure extends below axis line
                               ymax = qp$thresholds$prob, color = color))
        #Because threshline extended below axis limit, axis automatically scaled below it.
        #Restore it to its former values
    	p <- p + coord_cartesian(ylim=axisYrange)
      }
      if (ci) {
        qp$thresholdsci[[color]] <- factor(qp$thresholdsci[[color]])
        p <- p + geom_errorbarh(data = qp$thresholdsci,
                       height = .03, aes_string(x = 'threinf', xmin = 'threinf',
                       color = color, xmax = 'thresup', y = qp$thresholds$prob))
      }
    }
    else {
      if (averages) p <- p + geom_point(data = qp$averages,
                        aes_string(x = qp$x, y = 'prob'))
      if (curves) p <- p + geom_line(data = qp$curves,
                         aes_string(x = 'x', y = 'y'))
      if (thresholds) p <- p + geom_linerange(data = qp$thresholds,
                        aes_string(x = 'thre', ymin = qp$guess,
                               ymax = qp$thresholds$prob))
      if (ci) p <- p + geom_errorbarh(data = qp$thresholdsci,
                       height = .03, aes_string(x = 'threinf', xmin = 'threinf',
                       xmax = 'thresup', y = qp$thresholds$prob))
    }
  }
  p
}
