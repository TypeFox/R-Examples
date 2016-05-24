#' Plot the values of the parameters
#'
#' \code{plotpar_} is the standard evaluation SE function associated
#' to the non-standard evaluation NSE function \code{plotpar}.
#' \href{http://adv-r.had.co.nz/Computing-on-the-language.html}{SE functions can be more easily called from other functions.}
#' In SE functions, you need to quote the names of the variables.
#' @param qp output from quickpsy.
#' @param x Name of the variable to displayed in the x-axis.
#' @param panel Name of the variable to be split in panels.
#' @param xpanel Name of the variable to be split in horizontal panels.
#' @param ypanel Name of the variable to be split in vertical panels.
#' @param color Name of the variable codded by color.
#' @param geom If \code{'bar'} displays bars.
#' If \code{'point'} displays points (default is \code{'bar'}).
#' @param ci If \code{FALSE} confidence intervals are not plotted
#' (default is \code{TRUE}).
#' @seealso \code{\link{plotpar}}
#' @examples
#' library(MPDiR) # contains the Vernier data
#' fit <- quickpsy(Vernier, Phaseshift, NumUpward, N,
#'                 grouping = .(Direction, WaveForm, TempFreq), bootstrap = 'none')
#'
#' plotpar_(fit, x = 'WaveForm')
#' plotpar_(fit, xpanel = 'Direction')
#' plotpar_(fit, color = 'Direction')
#' plotpar_(fit, color = 'Direction', ypanel = 'WaveForm', geom = 'point')
#' @export
plotpar_ <- function(qp, x = NULL, panel = NULL, xpanel = NULL,
                           ypanel = NULL, color = NULL, geom = 'bar', ci  = T) {

  if (!('parci' %in% names(qp))) ci <- F

  p <- ggplot()
  if (qp$log) p <- p + ylab(paste0('log(', qp$x, ')'))
  else p <- p + ylab(qp$x)

  groups <- qp$groups
  ngroup <- length(groups)

  if (ngroup != 3) p <- p + facet_wrap(~parn, scales = 'free')

  if (ngroup == 0) { ###########################################################

    if (geom == 'bar') {
      p <- p + geom_bar(data = qp$par,
                        aes_string(x = 0, y = 'par'), fill = 'grey',
                        stat = 'identity', position = 'dodge') +
               theme(axis.title.x = element_blank(),
                        axis.text.x = element_blank())
      if (ci) p <- p + geom_errorbar(data = qp$parci, width = .5,
                 aes_string(x = 0, ymin = 'parinf', ymax = 'parsup'),
                 stat = 'identity', position = 'dodge')
    }
    if (geom == 'point') {
      p <- p + geom_point(data = qp$par,
                        aes_string(x = 0, y = 'par')) +
               theme(axis.title.x = element_blank(),
                        axis.text.x = element_blank())
      if (ci) p <- p + geom_linerange(data = qp$parci,
                       aes_string(x = 0, ymin = 'parinf', ymax = 'parsup'),
                       stat = 'identity', position = 'dodge')
    }
  }

  if (ngroup == 1) { ###########################################################
    if (is.null(color) && is.null(x)) color <- groups[[1]]

    if (!is.null(color)) {
      if (geom == 'bar') {
        qp$par[[color]] <- factor(qp$par[[color]])
        p <- p + geom_bar(data = qp$par,
                   aes_string(x = color,fill = color, y = 'par'),
                   stat = 'identity', position = 'dodge')
        if (ci) {
          qp$parci[[color]] <- factor(qp$parci[[color]])
          p <- p + geom_errorbar(data = qp$parci,
                   aes_string(x = color, ymin = 'parinf', ymax = 'parsup'),
                   stat = 'identity', position = 'dodge', width = .5)
        }
      }
      if (geom == 'point') {
        p <- p + geom_point(data = qp$par,
                          aes_string(x = color, color = color, y = 'par'))
        if (ci) p <- p + geom_linerange(data = qp$parci,
                   aes_string(x = color, color = color, ymin = 'parinf',
                   ymax = 'parsup'), stat = 'identity', position = 'dodge')
      }
    }
    if (is.null(color) && !is.null(x)) {
      if (geom == 'bar') {
        qp$par[[x]] <- factor(qp$par[[x]])
        p <- p + geom_bar(data = qp$par,
                          aes_string(x = x, y = 'par'), fill = 'grey',
                          stat = 'identity', position = 'dodge')
        if (ci) {
          qp$parci[[x]] <- factor(qp$parci[[x]])
          p <- p + geom_errorbar(data = qp$parci,
                         aes_string(x = x, width = .5, ymin = 'parinf',
                         ymax = 'parsup'), stat = 'identity',
                         position = 'dodge')
        }
      }
      if (geom == 'point') {
        p <- p + geom_point(data = qp$par,
                 aes_string(x = x, y = 'par')) +
                  geom_line(data = qp$par,
                                     aes_string(x = x, y = 'par'))
        if (ci) p <- p + geom_linerange(data = qp$parci,
                         aes_string(x = x, ymin = 'parinf', ymax = 'parsup'),
                         stat = 'identity', position = 'dodge')
      }
    }
  }

  if (ngroup == 2) { ###########################################################
    if (!is.null(x)) groups <- setdiff(groups, x)
    if (!is.null(color)) groups <- setdiff(groups, color)

    if (is.null(x)) {
      x <- groups[[1]]
      groups <- setdiff(groups, groups[[1]])
    }
    if (is.null(color)) color <- groups[[1]]
  }

  if (ngroup == 3) { ###########################################################
    if (!is.null(x)) groups <- setdiff(groups, x)
    if (!is.null(color)) groups <- setdiff(groups, color)
    if (!is.null(xpanel)) groups <- setdiff(groups, xpanel)
    if (!is.null(ypanel)) groups <- setdiff(groups, ypanel)
    if (!is.null(panel)) groups <- setdiff(groups, panel)

    if (is.null(x)) {
      x <- groups[[1]]
      groups <- setdiff(groups,groups[[1]])
    }
    if (is.null(color)) {
      color <- groups[[1]]
      groups <- setdiff(groups,groups[[1]])
    }
    if (is.null(xpanel) && is.null(ypanel) && is.null(panel)) {
      panel <- groups[[1]]
      p <- p + facet_grid(as.formula(paste0(panel, '~parn')),
                                   scales = 'free_y')
    }
    else {
      if (!is.null(xpanel)) p <- p +
        facet_grid(as.formula(paste0('parn~',xpanel)),
                            scales = 'free_y')

      if (!is.null(ypanel)) p <- p +
        facet_grid(as.formula(paste0(ypanel,'~parn')),
                            scales = 'free_y')
    }
  }

  if (ngroup == 2 || ngroup == 3) {
    qp$par[[color]] <- factor(qp$par[[color]])

    if (geom == 'bar') {
      qp$par[[x]] <- factor(qp$par[[x]])
      p <- p + geom_bar(data = qp$par,
               aes_string(x = x, fill = color, y = 'par'),
               stat = 'identity', position = 'dodge')
      if (ci) {
        qp$parci[[x]] <- factor(qp$parci[[x]])
        qp$parci[[color]] <- factor(qp$parci[[color]])
        p <- p + geom_errorbar(data = qp$parci, width = .5,
                 aes_string(x = x, fill = color, ymin = 'parinf',
                 ymax = 'parsup'), stat = 'identity',
                position = position_dodge(width=0.9))
      }
    }

    if (geom == 'point') {
      p <- p + geom_point(data = qp$par,
               aes_string(x = x, color = color, y = 'par')) +
               geom_line(data = qp$par, aes_string(x = x,
                        color = color, y = 'par', group =color))
      if (ci) {
        qp$parci[[color]] <- factor(qp$parci[[color]])
        p <- p + geom_linerange(data = qp$parci,
                               aes_string(x = x, color = color,
                                          ymin = 'parinf', ymax = 'parsup'))
      }
    }
  }
  p
}
