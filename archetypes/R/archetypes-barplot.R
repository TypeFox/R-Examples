

#' Barplot of archetypes.
#' @param height An \code{\link{archetypes}} object.
#' @param data The original data matrix.
#' @param which \code{below} creates a barplot for each archetype,
#'    \code{beside} creates one barplot with bars side by side.
#' @param which.beside Barplot according to \code{atypes} or \code{variables}.
#' @param which.below \code{compressed} plots the labels only once.
#' @param percentiles Show real values or percentile profiles.
#' @param below.compressed.height Height of additional tail subplot.
#' @param below.compressed.srt Rotations of the x-labels.
#' @param col.atypes Color of archetypes; only used in \code{below.compressed}.
#' @param ... Passed to the underlying \code{\link{barplot}} call.
#' @return Undefined.
#' @method barplot archetypes
#' @importFrom graphics barplot
#' @export
barplot.archetypes <- function(height, data,
                               which = c('below', 'beside'),
                               which.beside = c('atypes', 'variables'),
                               which.below = c('compressed', 'default'),
                               percentiles = FALSE,
                               below.compressed.height = 0.1,
                               below.compressed.srt = 0, 
                               col.atypes = NULL, ...) {

  ### Helpers:
  .beside.atypes <- function() {
    barplot(t(atypes), ylab=ylab, beside=TRUE, ylim=ylim, ...)
  }


  .beside.variables <- function() {
    barplot(atypes, ylab=ylab, beside=TRUE, ylim=ylim, ...)
  }


  .below.default <- function() {
    p <- nrow(atypes)

    layout(matrix(1:p, nrow = p, byrow = TRUE))
    for ( i in 1:p )
      barplot(atypes[i,], main=paste('Archetype', i),
              ylab=ylab, ylim=ylim, ...)
  }


  .below.compressed <- function() {
    p <- nrow(atypes) + 1
    heights <- c(rep(1, p - 1), below.compressed.height)

    layout(matrix(1:p, nrow = p, byrow = TRUE),
           heights = heights)
    for ( i in 1:(p - 1) ) {
      par(mar = c(0, 5, 1, 0) + 0.1)
      x.at <- barplot(atypes[i,], ylab = ylab, ylim = ylim,
                      names.arg = '', las = 2, col = col.atypes[i], ...)
      mtext(sprintf('Archetype %s', i), side = 2, line = 4,
            cex = par('cex'))
    }

    text(x.at, par("usr")[3] - 3, srt = below.compressed.srt,
         adj = 1, labels = colnames(atypes), xpd = NA)
  }


  .perc <- function(x, data, digits = 0) {
    Fn <- ecdf(data)
    round(Fn(x) * 100, digits = digits)
  }


  ### Plot:
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  which <- match.arg(which)

  if ( which == 'beside' )
    which.arg <- match.arg(which.beside)
  else
    which.arg <- match.arg(which.below)
   

  atypes <- parameters(height)
  rownames(atypes) <- sprintf('Archetype %s',
                              seq(length = nrow(atypes)))

  if ( !percentiles ) {
    ylab <- 'Value'
    ylim <- NULL
  }
  else {
    atypes <- sapply(seq(length = ncol(data)),
                     function(i)
                     .perc(atypes[, i], data[, i]))
    colnames(atypes) <- colnames(data)

    ylab <- 'Percentile'
    ylim <- c(0, 100)
  }

  do.call(sprintf('.%s.%s', which, which.arg), list())

  invisible(atypes)
}

