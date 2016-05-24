#' Plot weighted cumulative distribution functions
#'
#' Plot cumulative distribution functions, possibly broken down according to
#' conditioning variables and taking into account sample weights.
#'
#' Weights are directly extracted from the input object \code{inp} and are
#' taken into account by adjusting the step height.  To be precise, the
#' weighted step height for an observation is defined as its weight divided by
#' the sum of all weights\eqn{\ ( w_{i} / \sum_{j = 1}^{n} w_{j} ).}{.}
#'
#' @name spCdfplot
#' @aliases spCdfplot spCdfplot.default prepanelSpCdfplot panelSpCdfplot getCdf prepCdf prepCdf.data.frame prepCdf.default
#' @param inp an object of class \code{\linkS4class{simPopObj}} containing
#' survey sample and synthetic population data.
#' @param x a character vector specifying the columns of data available in the
#' sample and the population (specified in input object 'inp') to be plotted.
#' @param cond an optional character vector (of length 1, if used) specifying
#' the conditioning variable.
#' @param approx logicals indicating whether approximations of the cumulative
#' distribution functions should be computed.  The default is to use
#' \code{FALSE} for the survey data and \code{TRUE} for the population data.
#' @param n integers specifying the number of points at which the
#' approximations take place (see \code{\link[stats:approxfun]{approx}}).  It
#' is used wherever \code{approx} is \code{TRUE}.
#' @param bounds a logical indicating whether vertical lines should be drawn at
#' 0 and 1 (the bounds for cumulative distribution functions).
#' @param \dots further arguments to be passed to
#' \code{\link[lattice]{xyplot}}.
#' @return An object of class \code{"trellis"}, as returned by
#' \code{\link[lattice]{xyplot}}.
#' @author Andreas Alfons
#' @seealso \code{\link{spCdf}}, \code{\link[lattice]{xyplot}}
#' @keywords hplot
#' @export
#' @examples
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
#'   regModel = ~rb090+hsize+pl030+pb220a,
#'   upper=200000, equidist=FALSE, nr_cpus=1)
#' class(eusilcM)
#'
#' # plot results
#' spCdfplot(eusilcM, "netIncome", cond=NULL)
#' spCdfplot(eusilcM, "netIncome", cond="rb090", layout=c(1,2))
#' }
spCdfplot <- function(inp, x, cond = NULL, approx = c(FALSE, TRUE),
                      n = 10000, bounds = TRUE, ...) {
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

  # check 'approx'
  if(!is.logical(approx) || length(approx) == 0) approx <- formals()$approx
  else approx <- sapply(rep(approx, length.out=2), isTRUE)
  # check 'n'
  if(any(approx) && (!is.numeric(n) || length(n) == 0)) {
    stop("'n' is not numeric or does not have positive length")
  } else n <- ifelse(approx, n[1], NA)
  # check 'bounds'
  bounds <- isTRUE(bounds)
  # define labels for grouping variable
  lab <- c("Sample", "Population")

  ## construct objects for 'xyplot'
  # from sample
  tmp <- getCdf(x, weights.samp, cond, dataS, approx=approx[1], n=n[1], name=lab[1])
  values <- tmp$values
  app <- t(tmp$approx)
  # from population
  tmp <- getCdf(x, weights.pop, cond, dataP, approx=approx[2], n=n[2], name=lab[2])
  values <- rbind(values, tmp$values)
  app <- rbind(app, tmp$approx)
  ## construct formula for 'xyplot'
  form <- ".y~.x"  # basic formula
  if(length(x) > 1) cond <- c(".var", cond)
  if(!is.null(cond)) {
    cond <- paste(cond, collapse = " + ")  # conditioning variabels
    form <- paste(form, cond, sep=" | ")  # add conditioning to formula
  }
  ## define local version of 'xyplot'
  localXyplot <- function(form, values, xlab = NULL, ylab = NULL,
                          auto.key = TRUE, ...,
                          # these arguments are defined so that they aren't supplied twice:
                          x, data, allow.multiple, outer, panel, prepanel, groups) {
    # prepare legend
    if(isTRUE(auto.key)) auto.key <- list(points=FALSE, lines=TRUE)
    else if(is.list(auto.key)) {
      if(is.null(auto.key$points)) auto.key$points <- FALSE
      if(is.null(auto.key$lines)) auto.key$lines <- TRUE
    }
    command <- paste("xyplot(form, data=values, groups=.name,",
                     "panel=panelSpCdfplot, prepanel=prepanelSpCdfplot,",
                     "xlab=xlab, ylab=ylab, auto.key=auto.key, ...)")
    eval(parse(text=command))
  }
  ## call 'xyplot'
  localXyplot(as.formula(form), values, approx=app, bounds=bounds, ...)
}
NULL


## panel function
#' @rdname spCdfplot
#' @export
panelSpCdfplot <- function(x, y, approx, bounds = TRUE, ...) {
  if(isTRUE(bounds)) {
    panel.refline(h=0, ...)
    panel.refline(h=1, ...)
  }
  localPanelXyplot <- function(..., approx, type, distribute.type) {
    i <- packet.number()
    type <- ifelse(approx[,i], "l", "s")
    panel.xyplot(..., type=type, distribute.type=TRUE)
  }
  localPanelXyplot(x, y, approx=approx, ...)
}
NULL

## prepanel function
#' @rdname spCdfplot
#' @export
prepanelSpCdfplot <- function(x, y, ...) list(ylim=c(0,1))
NULL

## internal utility functions

# get data.frame and logical indicating approximation
#' @rdname spCdfplot
#' @export
#' @keywords internal
getCdf <- function(x, weights = NULL, cond = NULL, data, ..., name = "") {
  if ( is.null(cond) ) {
    x <- data[, x, with=FALSE]
    if ( length(weights) == 0 ) {
      w <- NULL
    } else {
      w <- data[[weights]]
    }
    prepCdf(x, w, ..., name=name)
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
      res <- prepCdf(x, w, ..., name=name)
      res$values <- cbind(res$values, rep(data[[cond]][1], nrow(res$values)))
      colnames(res$values)[ncol(res$values)] <- cond
      res
    })
    values <- do.call(rbind, lapply(tmp, function(x) x$values))
    approx <- as.vector(sapply(tmp, function(x) x$approx))
    list(values=values, approx=approx)
  }
}
NULL

# prepare one or more variables
#' @rdname spCdfplot
#' @export
#' @keywords internal
prepCdf <- function(x, w, ..., name = "") UseMethod("prepCdf")
NULL

#' @rdname spCdfplot
#' @export
#' @keywords internal
prepCdf.data.frame <- function(x, w, ..., name = "") {
  tmp <- lapply(x, prepCdf, w, ..., name=name)
  values <- mapply(function(x, v) cbind(x$values, .var=v), tmp, names(x), SIMPLIFY=FALSE, USE.NAMES=FALSE)
  values <- do.call(rbind, values)
  approx <- sapply(tmp, function(x) x$approx)
  list(values=values, approx=approx)
}
NULL

#' @rdname spCdfplot
#' @export
#' @keywords internal
prepCdf.default <- function(x, w, ..., name = "") {
  tmp <- spCdf(x, w, ...)
  values <- data.frame(.x=c(tmp$x[1], tmp$x), .y=c(0, tmp$y), .name=name)
  list(values=values, approx=tmp$approx)
}
NULL
