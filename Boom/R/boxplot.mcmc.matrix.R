BoxplotMcmcMatrix <- function(X, ylim = range(X), col.names,
                              row.names, truth, colors = NULL,
                              las = 0, ...) {
  ## Creates an organized set of boxplots for a sequence of MCMC draws
  ## of a matrix.
  ##
  ## Args:
  ##   X: A 3-way array representing the draws of the matrix to be
  ##     plotted.  The first index is the observation number.  The
  ##     second and third indices correspond to rows and columns of
  ##     the simulated matrices.
  ##   ylim:  Limits for the vertical axes.
  ##   col.names: Names for the column headings in the plot.  If
  ##     missing, these will be taken from the dimnames of X.
  ##   row.names: Names for the row headings in the plot.  If missing,
  ##     these will be taken from the dimnames of X.
  ##   truth: Either a scalar or a matrix with dimensions matching
  ##     X[1, , ], giving reference values to be plotted for each
  ##     coordinate.
  ##   colors:  A vector of colors to use for the boxes.
  ##   ...:  Extra arguments passed to boxplot.
  stopifnot(is.array(X) && length(dim(X)) == 3)
  nr <- dim(X)[2]
  nc <- dim(X)[3]
  have.truth <- FALSE
  if (!missing(truth)) {
    have.truth <- TRUE
    if (length(truth) == 1) {
      truth <- matrix(truth, nrow = nr, ncol = nc)
    }
    else if ((nrow(truth) != nr) || (ncol(truth) != nc)) {
      stop("dimension of  'truth' does not match dimension of 'X'\n")
    }
  }
  if (missing(row.names)) {
    row.names <- dimnames(X)[[2]]
  }
  if (is.null(row.names)) {
    row.names <- as.character(1:nr)
  }
  if (length(row.names) != nr) {
    stop("need ", nr, " entries in 'row.names'.  Have ", length(row.names), ".")
  }
  if (missing(col.names)) {
    col.names <- dimnames(X)[[3]]
  }
  if (is.null(col.names)) {
    col.names <- paste(1:nc)
  }
  if (length(col.names) != nc) {
    stop("need ", nc, " entries in 'col.names'.  Have ", length(col.names), ".")
  }
  opar <- par(mfrow = c(nr, 1), mar = rep(0, 4), oma = c(4, 8, 4, 4), las = las)
  on.exit(par(opar))
  rng <- range(X)
  for (j in 1:nr) {
    boxplot(split(X[, j, ], col(X[, j, ])), ylim = ylim,
            axes = FALSE, col = colors, ...)
    if (have.truth) {
      AddSegments(1:nc, truth[j, ], lwd = 3)
    }
    box()
    if (j %% 2) {
      axis(side = 2)
    } else {
      axis(side = 4)
    }
    mtext(row.names[j], side = 2, cex = par("cex.lab"),
          font = par("font.lab"), line = 2.5)
    if (j == nr) {
      axis(side = 1, at = 1:nc, labels = col.names)
    }
  }
}
