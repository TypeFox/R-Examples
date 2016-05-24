#' Weighted box plots
#'
#' Produce box-and-whisker plots of continuous or semi-continuous variables,
#' possibly broken down according to conditioning variables and taking into
#' account sample weights.
#'
#' Missing values are ignored for producing box plots and weights are directly
#' extracted from the input object \code{inp}.
#'
#' @name spBwplot
#' @aliases spBwplot panelSpBwplot getBwplotStats prepBwplotStats.data.frame prepBwplotStats.default
#' @param inp an object of class \code{\linkS4class{simPopObj}} containing
#' survey sample and synthetic population data.
#' @param x a character vector specifying the columns of data available in the
#' sample and the population (specified in input object 'inp') to be plotted.
#' @param cond an optional character vector (of length 1, if used) specifying
#' the conditioning variable.
#' @param horizontal a logical indicating whether the boxes should be
#' horizontal or vertical.
#' @param coef a numeric value that determines the extension of the whiskers.
#' @param zeros a logical indicating whether the variables specified by
#' \code{x} are semi-continuous, i.e., contain a considerable amount of zeros.
#' If \code{TRUE}, the box widths correspond to the proportion of non-zero data
#' points and the (weighted) box plot statistics are computed for these
#' non-zero data points only.
#' @param minRatio a numeric value in \eqn{(0,1]}; if \code{zeros} is
#' \code{TRUE}, the boxes may become unreadable for a large proportion of
#' zeros. In such a case, this can be used to specify a minimum ratio for the
#' box widths. Variable box widths for semi-continuous variables can be
#' suppressed by setting this value to 1.
#' @param do.out a logical indicating whether data points that lie beyond the
#' extremes of the whiskers should be plotted. Note that this is \code{FALSE}
#' by default.
#' @param \dots further arguments to be passed to
#' \code{\link[lattice:xyplot]{bwplot}}.
#' @return An object of class \code{"trellis"}, as returned by
#' \code{\link[lattice:xyplot]{bwplot}}.
#' @author Andreas Alfons and Bernhard Meindl
#' @seealso \code{\link{spBwplotStats}}, \code{\link[lattice:xyplot]{bwplot}}
#' @keywords hplot
#' @export
#' @examples
#'
#' ## these take some time and are not run automatically
#' ## copy & paste to the R command line
#'
#' set.seed(1234)  # for reproducibility
#' data(eusilcS)   # load sample data
#' \dontrun{
#' ## approx. 20 seconds computation time
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize",
#'   strata="db040", weight="db090")
#' simPop <- simStructure(data=inp, method="direct",
#'   basicHHvars=c("age", "rb090", "hsize", "pl030", "pb220a"))
#'
#' # multinomial model with random draws
#' eusilcM <- simContinuous(simPop, additional="netIncome",
#'   regModel  = ~rb090+hsize+pl030+pb220a+hsize,
#'   upper=200000, equidist=FALSE, nr_cpus=1)
#' class(eusilcM)
#'
#' # plot results
#' spBwplot(eusilcM, x="netIncome", cond=NULL)
#' spBwplot(eusilcM, x="netIncome", cond="rb090", layout=c(1,2))
#' }
#'
spBwplot <- function(inp, x, cond = NULL, horizontal = TRUE,
                     coef = 1.5, zeros = TRUE, minRatio = NULL,
                     do.out = FALSE, ...) {

  ## initializations
  if ( !class(inp) == "simPopObj" ) {
    stop("input argument 'inp' must be of class 'simPopObj'!\n")
  }

  weights.pop <- inp@pop@weight
  weights.samp <- inp@sample@weight
  dataS <- inp@sample@data
  dataP <- inp@pop@data

  if ( !is.character(x) || length(x) == 0 ) {
    stop("'x' must be a character vector of positive length!\n")
  }
  if ( !(all(x %in% colnames(dataP)) & (all(x %in% colnames(dataS)))) ) {
    stop("The variable names specified in argument 'x' must be available both in the population and the sample!\n")
  }

  if ( !is.null(cond) && !is.character(cond) ) {
    stop("'cond' must be a character vector or NULL!\n")
    if ( length(cond) != 1 ) {
      stop("argument 'cond' must have length 1!\n")
    }
  }
  if ( !(all(cond %in% colnames(dataP)) & (all(cond %in% colnames(dataS)))) ) {
    stop("The variable names specified in argument 'cond' must be available both in the population and the sample!")
  }

  horizontal <- isTRUE(horizontal)
  zeros <- isTRUE(zeros)
  do.out <- isTRUE(do.out)
  lab <- c("Sample", "Population")

  ## compute statistics for boxplots and construct objects for 'bwplot'
  # from sample
  tmp <- getBwplotStats(x, weights.samp, cond, dataS, coef=coef, zeros=zeros, do.out=do.out, name=lab[1])
  values <- tmp$values
  n <- t(tmp$n)
  nzero <- ifelse(zeros, t(tmp$nzero), NULL)
  out <- tmp$out
  # from population(s)
  tmp <- getBwplotStats(x, weights.pop, cond, dataP, coef=coef, zeros=zeros, do.out=do.out, name=lab[2])
  values <- rbind(values, tmp$values)
  n <- rbind(n, tmp$n)
  if ( zeros ) {
    nzero <- rbind(nzero, tmp$nzero)
  }
  out <- c(out, tmp$out)

  ## construct formula for 'bwplot'
  form <- ifelse(horizontal, ".name~.x", ".x~.name")  # basic formula
  if ( length(x) > 1 ) {
    cond <- c(".var", cond)
  }
  if ( !is.null(cond) ) {
    cond <- paste(cond, collapse = " + ")  # conditioning variabels
    form <- paste(form, cond, sep=" | ")  # add conditioning to formula
  }
  ## in case of semi-continuous variables define box widths
  if ( zeros ) {
    ratio <- n/(n+nzero)
    if ( !is.null(minRatio) ) {
      if ( !is.numeric(minRatio) || length(minRatio) != 1 || minRatio <= 0 || minRatio > 1 ) {
        stop("'minRatio' must be a single numeric value in [0,1]!\n")
      }
      ratio[ratio < minRatio] <- minRatio
    }
  } else {
    ratio <- NULL
  }
  ## define local version of 'bwplot'
  localBwplot <- function(form, values, xlab = NULL, ylab = NULL, ...,
                          # these arguments are defined so that they aren't supplied twice:
                          x, data, allow.multiple, outer, panel, groups) {
    bwplot(form, data=values, panel=panelSpBwplot, xlab=xlab, ylab=ylab, ...)
  }
  ## call 'bwplot'
  localBwplot(as.formula(form), values, horizontal=horizontal, coef=coef, zeros=zeros, ratio=ratio, do.out=FALSE, outliers=out, ...)
}
NULL


## panel function
#' @rdname spBwplot
#' @export
panelSpBwplot <- function(x, y, coef=1.5, zeros = TRUE, ratio, outliers, subscripts, ...) {
  out <- outliers[subscripts]
  if ( zeros ) {
    i <- packet.number()
    localPanelBwplot <- function(..., ratio, box.ratio, box.width) {
      panel.bwplot(..., box.ratio=ratio)
    }
    localPanelBwplot(x[!out], y[!out], coef=coef, ratio=ratio[,i], ...)
  } else {
    panel.bwplot(x[!out], y[!out], coef=coef, ...)
  }
  panel.points(x[out], y[out], ...)
}
NULL


## internal utility functions

# get data.frame and all required statistics
#' @rdname spBwplot
#' @export
#' @keywords internal
getBwplotStats <- function(x, weights = NULL, cond = NULL, data, ..., name = "") {
  if ( is.null(cond) ) {
    x <- data[, x, with=FALSE]
    if ( length(weights) == 0 ) {
      w <- NULL
    } else {
      w <- data[[weights]]
    }
    prepBwplotStats(x, w, ..., name=name)
  } else {
    spl <- split(data, data[[cond]])
    tmp <- lapply(spl, function(z) {
      data <- z
      x <- data[, x, with=FALSE]
      if ( length(weights) == 0 ) {
        w <- NULL
      } else {
        w <- data[[weights]]
      }
      res <- prepBwplotStats(x, w, ..., name=name)
      res$values <- cbind(res$values, rep(data[[cond]][1], nrow(res$values)))
      colnames(res$values)[ncol(res$values)] <- cond
      res
    })
    values <- do.call(rbind, lapply(tmp, function(x) x$values))
    n <- as.vector(sapply(tmp, function(x) x$n))
    nzero <- as.vector(sapply(tmp, function(x) x$nzero))
    out <- do.call(c, lapply(tmp, function(x) x$out))
    list(values=values, n=n, nzero=nzero, out=out)
  }
}
NULL

# prepare one or more variables
#' @rdname spBwplot
#' @export
#' @keywords internal
prepBwplotStats <- function(x, w, ..., name = "") UseMethod("prepBwplotStats")
NULL

#' @rdname spBwplot
#' @export
#' @keywords internal
prepBwplotStats.data.frame <- function(x, w, ..., name = "") {
  tmp <- lapply(x, prepBwplotStats, w, ..., name=name)
  values <- mapply(function(x, v) cbind(x$values, .var=v), tmp, names(x), SIMPLIFY=FALSE, USE.NAMES=FALSE)
  values <- do.call(rbind, values)
  n <- sapply(tmp, function(x) x$n)
  nzero <- sapply(tmp, function(x) x$nzero)
  out <- do.call(c, lapply(tmp, function(x) x$out))
  list(values=values, n=n, nzero=nzero, out=out)
}
NULL

#' @rdname spBwplot
#' @export
#' @keywords internal
prepBwplotStats.default <- function(x, w, ..., name = "") {
  stats <- spBwplotStats(x, w, ...)
  x <- c(stats$stats, stats$out)
  n <- stats$n
  nzero <- stats$nzero
  nstats <- length(stats$stats)
  nout <- length(stats$out)
  out <- rep.int(c(FALSE, TRUE), c(nstats, nout))
  values <- data.frame(.x=x, .name=name)
  list(values=values, n=n, nzero=nzero, out=out)
}
