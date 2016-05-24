# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' Plot a sequence of regression models
#'
#' Produce a plot of the coefficients, the values of the optimality criterion,
#' or diagnostic plots for a sequence of regression models, such as submodels
#' along a robust or groupwise least angle regression sequence, or sparse least
#' trimmed squares regression models for a grid of values for the penalty
#' parameter.
#'
#' @method plot seqModel
#' @aliases plot.rlars plot.grplars plot.tslarsP
#'
#' @param x  the model fit to be plotted.
#' @param p  an integer giving the lag length for which to produce the plot
#' (the default is to use the optimal lag length).
#' @param method  a character string specifying the type of plot.  Possible
#' values are \code{"coefficients"} to plot the coefficients from the submodels
#' via \code{\link{coefPlot}} (only for the \code{"seqModel"} and
#' \code{"sparseLTS"} methods), \code{"crit"} to plot the values of the
#' optimality criterion for the submodels via \code{\link{critPlot}}, or
#' \code{"diagnostic"} for diagnostic plots via \code{\link{diagnosticPlot}}.
#' @param \dots  additional arguments to be passed down.
#'
#' @return
#' An object of class \code{"ggplot"} (see \code{\link[ggplot2]{ggplot}}).
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{coefPlot}}, \code{\link{critPlot}},
#' \code{\link{diagnosticPlot}}, \code{\link{rlars}}, \code{\link{grplars}},
#' \code{\link{rgrplars}}, \code{\link{tslarsP}}, \code{\link{rtslarsP}},
#' \code{\link{tslars}}, \code{\link{rtslars}}, \code{\link{sparseLTS}}
#'
#' @example inst/doc/examples/example-plot.R
#'
#' @keywords hplot
#'
#' @export

plot.seqModel <- function(x, method = c("coefficients", "crit", "diagnostic"),
                          ...) {
  ## initializations
  method <- if(is.null(getSOpt(x))) "diagnostic" else match.arg(method)
  ## call plot function
  if(method == "coefficients") coefPlot(x, ...)
  else if(method == "crit") critPlot(x, ...)
  else diagnosticPlot(x, ...)
}


#' @rdname plot.seqModel
#' @method plot perrySeqModel
#' @export

plot.perrySeqModel <- function(x, method = c("crit", "diagnostic"), ...) {
  ## initializations
  method <- match.arg(method)
  ## call plot function
  if(method == "crit") critPlot(x, ...)
  else diagnosticPlot(x, ...)
}


#' @rdname plot.seqModel
#' @method plot tslars
#' @export

plot.tslars <- function(x, p, method = c("coefficients", "crit", "diagnostic"),
                        ...) {
  ## initializations
  method <- match.arg(method)
  ## call plot function
  if(method == "coefficients") coefPlot(x, p=p, ...)
  else if(method == "crit") critPlot(x, p=p, ...)
  else diagnosticPlot(x, p=p, ...)
}


#' @rdname plot.seqModel
#' @method plot sparseLTS
#' @export

plot.sparseLTS <- function(x, method = c("coefficients", "crit", "diagnostic"),
                           ...) {
  ## initializations
  crit <- x$crit
  if(is.null(crit) || inherits(crit, "fitSelect")) method <- "diagnostic"
  else method <- match.arg(method)
  ## call plot function
  if(method == "coefficients") coefPlot(x, ...)
  else if(method == "crit") critPlot(x, ...)
  else diagnosticPlot(x, ...)
}


#' @rdname plot.seqModel
#' @method plot perrySparseLTS
#' @export

plot.perrySparseLTS <- function(x, method = c("crit", "diagnostic"), ...) {
  ## initializations
  method <- match.arg(method)
  ## call plot function
  if(method == "crit") critPlot(x, ...)
  else diagnosticPlot(x, ...)
}

# ----------------------

## supplement the coefficients in a model with other useful information
## returns a data frame suitable for plotting with ggplot2

coefify <- function(model, ...) UseMethod("coefify")

coefify.seqModel <- function(model, zeros = FALSE, labels, ...) {
  # prepare coefficients and labels
  coef <- removeIntercept(t(coef(model, s=NULL)))
  sigmaX <- model$sigmaX
  if(!isTRUE(zeros)) {
    keep <- apply(coef != 0, 2, any)
    coef <- coef[, keep, drop=FALSE]
    sigmaX <- sigmaX[keep]
    if(!is.null(labels)) labels <- labels[keep]
  }
  # standardize coefficients
  coef <- sweep(coef, 2, sigmaX, "/", check.margin=FALSE)
  # prepare other information
  m <- ncol(coef)          # number of variables
  steps <- model$s         # step numbers
  nsteps <- length(steps)  # number of steps
  df <- model$df           # degrees of freedom
  vn <- colnames(coef)     # variable names
  # build data frame
  coefData <- data.frame(step=rep.int(steps, m),
                         df=rep.int(df, m), coefficient=as.numeric(coef),
                         variable=factor(rep(vn, each=nsteps), levels=vn))
  if(!is.null(labels))
    coefData$label <- rep(as.character(labels), each=nsteps)
  coefData
}

coefify.sparseLTS <- function(model, fit = c("reweighted", "raw", "both"),
                              zeros = FALSE, labels, ...) {
  # initializations
  fit <- match.arg(fit)
  zeros <- isTRUE(zeros)
  coef <- removeIntercept(t(coef(model, s=NULL, fit=fit)))
  df <- getComponent(model, "df", s=NULL, fit=fit)
  # prepare coefficients and labels
  if(!zeros) {
    keep <- apply(coef != 0, 2, any)
    coef <- coef[, keep, drop=FALSE]
    if(!is.null(labels)) labels <- labels[keep]
  }
  # check if predictor data is available to compute scale estimates
  if(is.null(x <- model$x)) {
    x <- try(model.matrix(model$terms), silent=TRUE)
    if(inherits(x, "try-error"))
      stop("scale estimates of predictor variables not available")
  }
  x <- removeIntercept(x)
  if(!zeros) x <- x[, keep, drop=FALSE]
#   # obtain scale estimates for predictors
#   lambda <- model$lambda      # tuning parameters
#   steps <- seq_along(lambda)  # step numbers
#   if(fit %in% c("reweighted", "both")) {
#     cdelta <- model$cnp2
#     wt <- as.matrix(wt(model, s=NULL, fit="raw"))
#     sigmaX <- do.call(rbind,
#                       lapply(steps, function(s) {
#                         xOk <- x[wt[, s] == 1, , drop=FALSE]
#                         apply(xOk, 2, sd) * cdelta[s]
#                       }))
#   } else sigmaX <- NULL
#   if(fit %in% c("raw", "both")) {
#     cdelta <- model$raw.cnp2
#     best <- as.matrix(model$best)
#     raw.sigmaX <- do.call(rbind,
#                           lapply(steps, function(s) {
#                             xBest <- x[best[, s], , drop=FALSE]
#                             apply(xBest, 2, sd) * cdelta
#                           }))
#   } else raw.sigmaX <- NULL
#   sigmaX <- rbind(sigmaX, raw.sigmaX)
#   # standardize coeffients
#   coef <- coef / sigmaX
#   # prepare other information
#   m <- ncol(coef)        # number of variables
#   sMax <- length(steps)  # number of steps
#   vn <- colnames(coef)   # variable names
#   obtain scale estimates for predictors
  n <- nrow(x)
  sigmaX <- apply(x, 2, function(x) {
    # standardize data
    xs <- robStandardize(x, fallback=TRUE)
    # detect good data points
    ok <- which(abs(xs) < qnorm(0.9875))
    nOk <- length(ok)
    # compute consistency factor
    if(nOk < n) {
      qn <- qnorm((nOk+n)/ (2*n))  # quantile for consistency factor
      cdelta <- 1 / sqrt(1-(2*n)/(nOk/qn)*dnorm(qn))
    } else cdelta <- 1  # consistency factor not necessary
    # compute standard deviation of good data points and multiply with
    # consistency factor
    sd(x[ok]) * cdelta
  })
  # standardize coeffients
  coef <- sweep(coef, 2, sigmaX, "/", check.margin=FALSE)
  # prepare other information
  m <- ncol(coef)             # number of variables
  lambda <- model$lambda      # tuning parameters
  steps <- seq_along(lambda)  # step numbers
  sMax <- length(steps)       # number of steps
  vn <- colnames(coef)        # variable names
  # build data frame
  if(fit == "both") {
    fits <- c("reweighted", "raw")
    coefData <- data.frame(
      fit=rep.int(factor(rep(fits, each=sMax), levels=fits), m),
      lambda=rep.int(lambda, 2*m), step=rep.int(steps, 2*m),
      df=rep.int(df, m), coefficient=as.numeric(coef),
      variable=factor(rep(vn, each=2*sMax), levels=vn))
    if(!is.null(labels))
      coefData$label <- rep(as.character(labels), each=2*sMax)
  } else {
    coefData <- data.frame(
      lambda=rep.int(lambda, m), step=rep.int(steps, m),
      df=rep.int(df, m), coefficient=as.numeric(coef),
      variable=factor(rep(vn, each=sMax), levels=vn))
    if(!is.null(labels))
      coefData$label <- rep(as.character(labels), each=sMax)
  }
  coefData
}


#' Coefficient plot of a sequence of regression models
#'
#' Produce a plot of the coefficients from a sequence of regression models,
#' such as submodels along a robust or groupwise least angle regression
#' sequence, or sparse least trimmed squares regression models for a grid of
#' values for the penalty parameter.
#'
#' @aliases coefPlot.rlars coefPlot.grplars coefPlot.tslarsP
#'
#' @param x  the model fit to be plotted.
#' @param p  an integer giving the lag length for which to produce the plot
#' (the default is to use the optimal lag length).
#' @param fit  a character string specifying for which estimator to produce the
#' plot.  Possible values are \code{"reweighted"} (the default) for the
#' reweighted fits, \code{"raw"} for the raw fits, or \code{"both"} for both
#' estimators.
#' @param abscissa  a character string specifying what to plot on the
#' \eqn{x}-axis.  Possible values are \code{"step"} for the step number (the
#' default), or \code{"df"} for the degrees of freedom.
#' @param zeros  a logical indicating whether predictors that never enter the
#' model and thus have zero coefficients should be included in the plot
#' (\code{TRUE}) or omitted (\code{FALSE}, the default).  This is useful if the
#' number of predictors is much larger than the number of observations, in
#' which case many coefficients are never nonzero.
#' @param size  a numeric vector of length three giving the line width, the
#' point size and the label size, respectively.
#' @param labels  an optional character vector containing labels for the
#' predictors.  Plotting labels can be suppressed by setting this to
#' \code{NULL}.
#' @param offset   an integer giving the offset of the labels from the
#' corresponding coefficient values from the last step (i.e., the number of
#' blank characters to be prepended to the label).
#' @param \dots  for the generic function, additional arguments to be passed
#' down to methods.  For the \code{"tslars"} method, additional arguments to be
#' passed down to the \code{"seqModel"} method.  For the other methods,
#' additional arguments to be passed down to \code{\link[ggplot2]{geom_line}}
#' and \code{\link[ggplot2]{geom_point}}.
#'
#' @return
#' An object of class \code{"ggplot"} (see \code{\link[ggplot2]{ggplot}}).
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link{rlars}},
#' \code{\link{grplars}}, \code{\link{rgrplars}}, \code{\link{tslarsP}},
#' \code{\link{rtslarsP}}, \code{\link{tslars}}, \code{\link{rtslars}},
#' \code{\link{sparseLTS}}
#'
#' @example inst/doc/examples/example-coefPlot.R
#'
#' @keywords hplot
#'
#' @export
#' @importFrom utils head tail

coefPlot <- function(x, ...) UseMethod("coefPlot")


#' @rdname coefPlot
#' @method coefPlot seqModel
#' @export

coefPlot.seqModel <- function(x, abscissa = c("step", "df"), zeros = FALSE,
                              size = c(0.5, 2, 4), labels, offset = 1, ...) {
  ## initializations
  if(missing(labels)) labels <- defaultLabels(x)  # default labels
  ## extract coefficient data extended with other information
  coefData <- coefify(x, zeros=zeros, labels=labels)
  ## construct data frame for labels
  maxStep <- max(coefData$step)
  labelData <- coefData[coefData$step == maxStep, ]
  ## call workhorse function
  ggCoefPlot(coefData, labelData, abscissa=abscissa, size=size,
             offset=offset, ...)
}


#' @rdname coefPlot
#' @method coefPlot tslars
#' @export

coefPlot.tslars <- function(x, p, ...) {
  ## check lag length
  if(missing(p) || !is.numeric(p) || length(p) == 0) p <- x$pOpt
  if(length(p) > 1) {
    warning("multiple lag lengths not yet supported")
    p <- p[1]
  }
  pMax <- x$pMax
  if(p < 1) {
    p <- 1
    warning("lag length too small, using lag length 1")
  } else if(p > pMax) {
    p <- pMax
    warning(sprintf("lag length too large, using maximum lag length %d", p))
  }
  ## call plot function for specified lag length
  coefPlot(x$pFit[[p]], ...)
}


#' @rdname coefPlot
#' @method coefPlot sparseLTS
#' @export

coefPlot.sparseLTS <- function(x, fit = c("reweighted", "raw", "both"),
                               abscissa = c("step", "df"), zeros = FALSE,
                               size = c(0.5, 2, 4), labels, offset = 1, ...) {
  ## initializations
  fit <- match.arg(fit)
  abscissa <- match.arg(abscissa)
  if(missing(labels)) labels <- defaultLabels(x)  # default labels
  ## extract coefficient data extended with other information
  coefData <- coefify(x, fit=fit, zeros=zeros, labels=labels)
  ## construct data frame for labels
  maxX <- max(coefData[, abscissa])
  labelData <- coefData[coefData[, abscissa] == maxX, ]
  if(abscissa == "df") {
    # maximum degree of freedom may occur in more than one step
    # ensure that label is only drawn once for largest step number
    by <- if(fit == "both") c("fit", "variable") else "variable"
    keep <- split(rownames(labelData), labelData[, by])
    keep <- sapply(keep, tail, 1)
    labelData <- labelData[keep, ]
  }
  ## call workhorse function
  p <- ggCoefPlot(coefData, labelData, abscissa=abscissa, size=size,
                  offset=offset, ...)
  if(fit == "both") {
    # split plot into different panels
    p <- p + facet_grid(. ~ fit)
  }
  p
}


## workhorse function
ggCoefPlot <- function(coefData, labelData, abscissa = c("step", "df"),
                       zeros = FALSE, size = c(0.5, 2, 4), labels, offset = 1,
                       main = NULL, xlab, ylab, ..., mapping, data) {
  # initializations
  abscissa <- match.arg(abscissa)
  size <- as.numeric(size)
  size <- c(size, rep.int(NA, max(0, 3-length(size))))[1:3]  # ensure length 3
  size <- ifelse(is.na(size), eval(formals()$size), size)    # fill NA's
  # define default axis labels
  if(missing(xlab))
    xlab <- switch(abscissa, step="Step", df="Degrees of freedom")
  if(missing(ylab)) ylab <- "Standardized coefficients"
  # define aesthetic mapping for plotting coefficients
  coefMapping <- aes_string(x=abscissa, y="coefficient", color="variable")
  # define aesthetic mapping for plotting x-axis grid and labels
  offset <- paste(rep.int(" ", offset), collapse="")  # whitespace
  labelData$label <- paste(offset, labelData$label, sep="")
  labelMapping <- aes_string(x=abscissa, y="coefficient", label="label")
  # draw minor grid lines for each step, but leave
  # major grid lines and tick marks pretty
  gridX <- unique(coefData[, abscissa])
  # create plot
  ggplot(coefData) +
    geom_line(coefMapping, size=size[1], ...) +
    geom_point(coefMapping, size=size[2], ...) +
    geom_text(labelMapping, data=labelData,
              hjust=0, size=size[3], alpha=0.4) +
    scale_x_continuous(minor_breaks=gridX) +
    theme(legend.position="none") +
    labs(title=main, x=xlab, y=ylab)
}

# ----------------------

#' Optimality criterion plot of a sequence of regression models
#'
#' Produce a plot of the values of the optimality criterion for a sequence of
#' regression models, such as submodels along a robust or groupwise least angle
#' regression sequence, or sparse least trimmed squares regression models for
#' a grid of values for the penalty parameter.
#'
#' @aliases critPlot.rlars critPlot.grplars critPlot.tslarsP
#'
#' @param x  the model fit to be plotted.
#' @param p  an integer giving the lag length for which to produce the plot
#' (the default is to use the optimal lag length).
#' @param fit  a character string specifying for which estimator to produce the
#' plot.  Possible values are \code{"reweighted"} (the default) for the
#' reweighted fits, \code{"raw"} for the raw fits, or \code{"both"} for both
#' estimators.
#' @param size  a numeric vector of length two giving the line width and the
#' point size, respectively.
#' @param \dots  for the generic function, additional arguments to be passed
#' down to methods.  For the \code{"tslars"} method, additional arguments
#' to be passed down to the \code{"seqModel"} method.  For the
#' \code{"seqModel"} and \code{"sparseLTS"} methods, additional arguments
#' to be passed down to \code{\link[ggplot2]{geom_line}} and
#' \code{\link[ggplot2]{geom_point}}.    For the \code{"perrySeqModel"} and
#' \code{"perrySparseLTS"} methods, additional arguments to be passed down
#' to the \code{\link[perry:perryPlot]{plot}} method for the prediction error
#' results.
#'
#' @return
#' An object of class \code{"ggplot"} (see \code{\link[ggplot2]{ggplot}}).
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[perry]{perryPlot}},
#' \code{\link{rlars}}, \code{\link{grplars}}, \code{\link{rgrplars}},
#' \code{\link{tslarsP}}, \code{\link{rtslarsP}}, \code{\link{tslars}},
#' \code{\link{rtslars}}, \code{\link{sparseLTS}}
#'
#' @example inst/doc/examples/example-critPlot.R
#'
#' @keywords hplot
#'
#' @export

critPlot <- function(x, ...) UseMethod("critPlot")


#' @rdname critPlot
#' @method critPlot seqModel
#' @export

critPlot.seqModel <- function(x, size = c(0.5, 2), ...) {
  ## extract information from object
  crit <- x$crit
  if(is.null(crit)) stop("optimality criterion data not available")
  ## construct data frame for ggplot2 graphics
  critData <- fortify(crit, data=data.frame(step=x$s))
  ## call workhorse function
  ggCritPlot(critData, abscissa="step", size=size, ...)
}


#' @rdname critPlot
#' @method critPlot perrySeqModel
#' @export

critPlot.perrySeqModel <- function(x, ...) {
  ## local plot function for prediction error results to override defaults
  localPlot <- function(x, method = c("line", "dot", "box", "density"),
                        xlab = "Step", ...) {
    # initializations
    if(x$splits$R == 1) {
      choices <- eval(formals()[["method"]])
      if(identical(method, choices)) method <- "line"
      else method <- match.arg(method, c("line", "dot"))
    } else method <- match.arg(method)
    # call perryPlot() for prediction error results
    p <- perryPlot(x, method=method, xlab=xlab, ...)
    if(method != "density") {
      p <- p + scale_x_continuous(minor_breaks=fits(x))
    }
    p
  }
  ## call local plot function
  localPlot(x, ...)
}


#' @rdname critPlot
#' @method critPlot tslars
#' @export

critPlot.tslars <- function(x, p, ...) {
  ## check lag length
  if(missing(p) || !is.numeric(p) || length(p) == 0) p <- x$pOpt
  if(length(p) > 1) {
    warning("multiple lag lengths not yet supported")
    p <- p[1]
  }
  pMax <- x$pMax
  if(p < 1) {
    p <- 1
    warning("lag length too small, using lag length 1")
  } else if(p > pMax) {
    p <- pMax
    warning(sprintf("lag length too large, using maximum lag length %d", p))
  }
  ## call plot function for specified lag length
  critPlot(x$pFit[[p]], ...)
}


#' @rdname critPlot
#' @method critPlot sparseLTS
#' @export

critPlot.sparseLTS <- function(x, fit = c("reweighted", "raw", "both"),
                               size = c(0.5, 2), ...) {
  ## initializations
  crit <- x$crit
  if(is.null(crit)) stop("optimality criterion data not available")
  fit <- match.arg(fit)
  select <- if(fit == "both") NULL else fit
  ## construct data frame for ggplot2 graphics
  critData <- fortify(crit, data=data.frame(lambda=x$lambda), select=select)
  ## call workhorse function
  p <- ggCritPlot(critData, abscissa="lambda", size=size, ...)
  if(fit == "both") {
    # split plot into different panels
    p <- p + facet_grid(. ~ fit)
  }
  p
}


#' @rdname critPlot
#' @method critPlot perrySparseLTS
#' @export

critPlot.perrySparseLTS <- function(x, fit = c("reweighted", "raw", "both"),
                                    ...) {
  ## local plot function for prediction error results to override defaults
  localPlot <- function(x, method = c("line", "dot", "box", "density"),
                        fit = select, select = "reweighted", xlab = "lambda",
                        ...) {
    # initializations
    if(x$splits$R == 1) {
      choices <- eval(formals()[["method"]])
      if(identical(method, choices)) method <- "line"
      else method <- match.arg(method, c("line", "dot"))
    } else method <- match.arg(method)
    # call perryPlot() for prediction error results
    if(is.null(fit)) p <- perryPlot(x, method=method, xlab=xlab, ...)
    else {
      p <- perryPlot(x, method=method, select=fit, facets=NULL, xlab=xlab, ...)
    }
    if(method != "density") {
      p <- p + scale_x_reverse(minor_breaks=x$tuning[, "lambda"])
    }
    p
  }
  ## call local plot function
  if(missing(fit)) localPlot(x, ...)
  else {
    fit <- match.arg(fit)
    if(fit == "both") fit <- NULL
    localPlot(x, fit=fit, ...)
  }
}


## workhorse function
ggCritPlot <- function(data, abscissa = c("index", "step", "lambda"),
                       size = c(0.5, 2), main = NULL, xlab, ylab, ...,
                       mapping) {
  # initializations
  abscissa <- match.arg(abscissa)
  crit <- setdiff(names(data), c("fit", "index", "step", "lambda"))
  size <- as.numeric(size)
  size <- c(size, rep.int(NA, max(0, 2-length(size))))[1:2]  # ensure length 2
  size <- ifelse(is.na(size), eval(formals()$size), size)    # fill NA's
  # define default axis labels
  if(missing(xlab)) {
    xlab <- switch(abscissa, index="Index", step="Step", lambda="lambda")
  }
  if(missing(ylab)) ylab <- crit
  # define aesthetic mapping for plotting coefficients
  mapping <- aes_string(x=abscissa, y=crit)
  # draw minor grid lines for each step, but leave
  # major grid lines and tick marks pretty
  gridX <- unique(data[, abscissa])
  # create plot
  scale_x <- if(abscissa == "lambda") scale_x_reverse else scale_x_continuous
  ggplot(data, mapping) +
    geom_line(size=size[1], ...) +
    geom_point(size=size[2], ...) +
    scale_x(minor_breaks=gridX) +
    labs(title=main, x=xlab, y=ylab)
}

# ----------------------

## construct data frame for labels based on some order
labelify <- function(data, which, id.n = NULL) {
  # initializations
  if(isTRUE(id.n < 1)) return(NULL)
  by <- intersect(c("step", "fit"), names(data))
  ord <- data[, which]
  if(which == "residual") ord <- abs(ord)
  if(is.null(id.n)) {
    # define outlier indicator
    out <- data[, "weight"] == 0                       # regression outliers
    if(which == "xyd") out <- out | data[, "leverage"] # all outliers
  }
  # find the id.n largest observations to keep
  # if NULL, id.n is computed as the number of outliers
  if(length(by) == 0) {
    # use the whole data set
    if(is.null(id.n)) id.n <- sum(out)
    if(id.n == 0) return(NULL)
    keep <- head(order(ord, decreasing=TRUE), id.n)
  } else {
    # split the data set according to the selected variables
    keep <- tapply(seq_len(nrow(data)), data[, by, drop=FALSE],
                   function(i) {
                     if(is.null(id.n)) id.n <- sum(out[i])
                     if(id.n > 0) {
                       largest <- head(order(ord[i], decreasing=TRUE), id.n)
                       i[largest]
                     }
                   })
    # combine indices to keep
    keep <- unlist(keep, use.names=FALSE)
    if(length(keep) == 0) return(NULL)
  }
  # return data frame with selected observations
  data[keep, ]
}


#' Diagnostic plots for a sequence of regression models
#'
#' Produce diagnostic plots for a sequence of regression models, such as
#' submodels along a robust least angle regression sequence, or sparse least
#' trimmed squares regression models for a grid of values for the penalty
#' parameter.  Four plots are currently implemented.
#'
#' In the normal Q-Q plot of the standardized residuals, a reference line is
#' drawn through the first and third quartile.  The \code{id.n} observations
#' with the largest distances from that line are identified by a label (the
#' observation number).  The default for \code{id.n} is the number of
#' regression outliers, i.e., the number of observations whose residuals are
#' too large (cf. \code{\link{wt}}).
#'
#' In the plots of the standardized residuals versus their index or the fitted
#' values, horizontal reference lines are drawn at 0 and +/-2.5.  The
#' \code{id.n} observations with the largest absolute values of the
#' standardized residuals are identified by a label (the observation
#' number).  The default for \code{id.n} is the number of regression outliers,
#' i.e., the number of observations whose absolute residuals are too large (cf.
#' \code{\link{wt}}).
#'
#' For the regression diagnostic plot, the robust Mahalanobis distances of the
#' predictor variables are computed via the MCD based on only those predictors
#' with non-zero coefficients (see
#' \code{\link[robustbase]{covMcd}}).  Horizontal reference lines are drawn at
#' +/-2.5 and a vertical reference line is drawn at the upper 97.5\% quantile
#' of the \eqn{\chi^{2}}{chi-squared} distribution with \eqn{p} degrees of
#' freedom, where \eqn{p} denotes the number of predictors with non-zero
#' coefficients.  The \code{id.n} observations with the largest absolute values
#' of the standardized residuals and/or largest robust Mahalanobis distances
#' are identified by a label (the observation number).  The default for
#' \code{id.n} is the number of all outliers: regression outliers (i.e.,
#' observations whose absolute residuals are too large, cf. \code{\link{wt}})
#' and leverage points (i.e., observations with robust Mahalanobis distance
#' larger than the 97.5\% quantile of the \eqn{\chi^{2}}{chi-squared}
#' distribution with \eqn{p} degrees of freedom).
#'
#' @aliases diagnosticPlot.rlars diagnosticPlot.grplars diagnosticPlot.tslarsP
#'
#' @param x  the model fit for which to produce diagnostic plots, or a data
#' frame containing all necessary information for plotting (as generated by the
#' corresponding \code{\link[=fortify.seqModel]{fortify}} method).
#' @param p  an integer giving the lag length for which to produce the plot
#' (the default is to use the optimal lag length).
#' @param s  for the \code{"seqModel"} method, an integer vector giving
#' the steps of the submodels  for which to produce diagnostic plots (the
#' default is to use the optimal submodel).  For the \code{"sparseLTS"} method,
#' an integer vector giving the indices of the models for which to produce
#' diagnostic plots (the default is to use the optimal model for each of the
#' requested fits).
#' @param fit  a character string specifying for which fit to produce
#' diagnostic plots.  Possible values are \code{"reweighted"} (the default) for
#' diagnostic plots for the reweighted fit, \code{"raw"} for diagnostic plots
#' for the raw fit, or \code{"both"} for diagnostic plots for both fits.
#' @param covArgs  a list of arguments to be passed to
#' \code{\link[robustbase]{covMcd}} for the regression diagnostic plot (see
#' \dQuote{Details}).
#' @param which  a character string indicating which plot to show.  Possible
#' values are \code{"all"} (the default) for all of the following, \code{"rqq"}
#' for a normal Q-Q plot of the standardized residuals, \code{"rindex"} for a
#' plot of the standardized residuals versus their index, \code{"rfit"} for a
#' plot of the standardized residuals versus the fitted values, or
#' \code{"rdiag"} for a regression diagnostic plot  (standardized residuals
#' versus robust Mahalanobis distances of the predictor variables).
#' @param ask  a logical indicating whether the user should be asked before
#' each plot (see \code{\link[grDevices]{devAskNewPage}}). The default is to
#' ask if all plots are requested and not ask otherwise.
#' @param facets  a faceting formula to override the default behavior.  If
#' supplied, \code{\link[ggplot2]{facet_wrap}} or
#' \code{\link[ggplot2]{facet_grid}} is called depending on whether the formula
#' is one-sided or two-sided.
#' @param size  a numeric vector of length two giving the point and label size,
#' respectively.
#' @param id.n  an integer giving the number of the most extreme observations
#' to be identified by a label.  The default is to use the number of identified
#' outliers, which can be different for the different plots.  See
#' \dQuote{Details} for more information.
#' @param \dots  for the generic function, additional arguments to be passed
#' down to methods.  For the \code{"tslars"} method, additional arguments to be
#' passed down to the \code{"seqModel"} method.  For the \code{"perrySeqModel"}
#' and \code{"perrySparseLTS"} method, additional arguments to be passed down
#' to the \code{"seqModel"} and \code{"sparseLTS"} method, respectively.  For
#' the \code{"seqModel"} and \code{"sparseLTS"} methods, additional arguments
#' to be passed down to the default method.  For the default method, additional
#' arguments to be passed down to \code{\link[ggplot2]{geom_point}}.
#'
#' @return
#' If only one plot is requested, an object of class \code{"ggplot"} (see
#' \code{\link[ggplot2]{ggplot}}), otherwise a list of such objects.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link{rlars}},
#' \code{\link{grplars}}, \code{\link{rgrplars}}, \code{\link{tslarsP}},
#' \code{\link{rtslarsP}}, \code{\link{tslars}}, \code{\link{rtslars}},
#' \code{\link{sparseLTS}}, \code{\link[robustbase:ltsPlot]{plot.lts}}
#'
#' @example inst/doc/examples/example-diagnosticPlot.R
#'
#' @keywords hplot
#'
#' @export
#' @import robustbase
#' @importFrom grDevices devAskNewPage

diagnosticPlot <- function(x, ...) UseMethod("diagnosticPlot")


#' @rdname diagnosticPlot
#' @method diagnosticPlot seqModel
#' @export

diagnosticPlot.seqModel <- function(x, s = NA, covArgs = list(), ...) {
  # call default method with all information required for plotting
  diagnosticPlot(fortify(x, s=s, covArgs=covArgs), ...)
}


#' @rdname diagnosticPlot
#' @method diagnosticPlot perrySeqModel
#' @export

diagnosticPlot.perrySeqModel <- function(x, ...) {
  # call method for component 'finalModel'
  diagnosticPlot(x$finalModel, ...)
}


#' @rdname diagnosticPlot
#' @method diagnosticPlot tslars
#' @export

diagnosticPlot.tslars <- function(x, p, ...) {
  ## check lag length
  if(missing(p) || !is.numeric(p) || length(p) == 0) p <- x$pOpt
  if(length(p) > 1) {
    warning("multiple lag lengths not yet supported")
    p <- p[1]
  }
  pMax <- x$pMax
  if(p < 1) {
    p <- 1
    warning("lag length too small, using lag length 1")
  } else if(p > pMax) {
    p <- pMax
    warning(sprintf("lag length too large, using maximum lag length %d", p))
  }
  ## call plot function for specified lag length
  diagnosticPlot(x$pFit[[p]], ...)
}


#' @rdname diagnosticPlot
#' @method diagnosticPlot sparseLTS
#' @export

diagnosticPlot.sparseLTS <- function(x, s = NA,
                                     fit = c("reweighted", "raw", "both"),
                                     covArgs = list(), ...) {
  # call default method with all information required for plotting
  diagnosticPlot(fortify(x, s=s, fit=fit, covArgs=covArgs), ...)
}


#' @rdname diagnosticPlot
#' @method diagnosticPlot perrySparseLTS
#' @export

diagnosticPlot.perrySparseLTS <- function(x, ...) {
  # call method for component 'finalModel'
  diagnosticPlot(x$finalModel, ...)
}


#' @rdname diagnosticPlot
#' @method diagnosticPlot default
#' @export

diagnosticPlot.default <- function(x, which = c("all", "rqq","rindex",
                                                "rfit", "rdiag"),
                                   ask = (which == "all"),
                                   facets = attr(x, "facets"),
                                   size = c(2, 4), id.n = NULL, ...) {
  # initializations
  which <- match.arg(which)
  size <- as.numeric(size)
  size <- c(size, rep.int(NA, max(0, 2-length(size))))[1:2]  # ensure length 2
  size <- ifelse(is.na(size), eval(formals()$size), size)    # fill NA's
  # call functions for selected plots
  if(which == "all") {
    oldAsk <- devAskNewPage(ask)  # ask for new page (if requested)
    on.exit(devAskNewPage(oldAsk))
    # residual Q-Q plot
    p <- try(rqqPlot(x, facets=facets, size=size, id.n=id.n, ...),
             silent=TRUE)
    if(inherits(p, "try-error")) {
      warn <- gsub("Error in", "In", p)
      warning(warn, call.=FALSE)
      res <- list()
    } else {
      print(p)
      res <- list(rqq=p)
    }
    # residuals vs indices plot
    p <- try(residualPlot(x, abscissa="index", facets=facets,
                          size=size, id.n=id.n, ...), silent=TRUE)
    if(inherits(p, "try-error")) {
      warn <- gsub("Error in", "In", p)
      warning(warn, call.=FALSE)
    } else {
      print(p)
      res$rindex <- p
    }
    # residuals vs fitted plot
    p <- try(residualPlot(x, abscissa="fitted", facets=facets,
                          size=size, id.n=id.n, ...), silent=TRUE)
    if(inherits(p, "try-error")) {
      warn <- gsub("Error in", "In", p)
      warning(warn, call.=FALSE)
    } else {
      print(p)
      res$rfit <- p
    }
    # regression diagnostic plot
    p <- try(rdiagPlot(x, facets=facets, size=size, id.n=id.n, ...),
             silent=TRUE)
    if(inherits(p, "try-error")) {
      warn <- gsub("Error in", "In", p)
      warning(warn, call.=FALSE)
    } else {
      print(p)
      res$rdiag <- p
    }
    invisible(res)
  } else if(which == "rqq") {
    # residual Q-Q plot
    rqqPlot(x, facets=facets, size=size, id.n=id.n, ...)
  } else if(which == "rindex") {
    # residuals vs indices plot
    residualPlot(x, abscissa="index", facets=facets,
                 size=size, id.n=id.n, ...)
  } else if(which == "rfit") {
    # residuals vs fitted plot
    residualPlot(x, abscissa="fitted", facets=facets,
                 size=size, id.n=id.n, ...)
  } else if(which == "rdiag") {
    # regression diagnostic plot
    rdiagPlot(x, facets=facets, size=size, id.n=id.n, ...)
  }
}

# ----------------------

rqqPlot <- function(data, facets = attr(data, "facets"), size = c(2, 4),
                    id.n = NULL, main, xlab, ylab, ..., mapping) {
  # define aesthetic mapping for Q-Q plot
  mapping <- aes_string(x="theoretical", y="residual", color="classification")
  # extract data frame for reference line
  lineData <- attr(data, "qqLine")
  # construct data frame for labels
  labelData <- labelify(data, which="qqd", id.n=id.n)
  # define default title and axis labels
  if(missing(main)) main <- "Normal Q-Q plot"
  if(missing(xlab)) xlab <- "Quantiles of the standard normal distribution"
  if(missing(ylab)) ylab <- "Standardized residual"
  # create plot
  p <- ggplot(data)
  if(!is.null(lineData)) {
    # add reference line
    lineMapping <- aes_string(intercept="intercept", slope="slope")
    p <- p + geom_abline(lineMapping, lineData, alpha=0.4)
  }
  p <- p + geom_point(mapping, size=size[1], ...)
  if(!is.null(labelData)) {
    # add labels for observations with largest distances
    labelMapping <- aes_string(x="theoretical", y="residual", label="index")
    p <- p + geom_text(labelMapping, data=labelData,
                       hjust=0, size=size[2], alpha=0.4)
  }
  p <- p + labs(title=main, x=xlab, y=ylab)
  if(!is.null(facets)) {
    # split plot into different panels
    if(length(facets) == 2) p <- p + facet_wrap(facets)
    else p <- p + facet_grid(facets)
  }
  p
}

## compute theoretical quantiles
qqNorm <- function(y) {
  # TODO: NA handling
  n <- length(y)                # number of observations
  prob <- ppoints(n)            # probabilities
  qnorm(prob)[order(order(y))]  # theoretical quantiles in original order
}

## compute intercept and slope of reference line
qqLine <- function(y) {
  prob <- c(0.25, 0.75)
  ly <- quantile(y, prob, na.rm=TRUE, names=FALSE)
  lx <- qnorm(prob)
  slope <- diff(ly) / diff(lx)
  intercept <- ly[1] - slope * lx[1]
  list(intercept=intercept, slope=slope)
}

# ----------------------

## plot standardized residuals vs indices or fitted values

residualPlot <- function(data, abscissa = c("index", "fitted"),
                         facets = attr(data, "facets"), size = c(2, 4),
                         id.n = NULL, main, xlab, ylab, ..., mapping) {
  ## initializations
  abscissa <- match.arg(abscissa)
  # define aesthetic mapping for residual plot
  mapping <- aes_string(x=abscissa, y="residual", color="classification")
  ## construct data frame for labels
  labelData <- labelify(data, which="residual", id.n=id.n)
  # define default title and axis labels
  if(missing(main)) {
    postfix <- switch(abscissa, index="indices", fitted="fitted values")
    main <- paste("Residuals vs", postfix)
  }
  if(missing(xlab)) {
    xlab <- switch(abscissa, index="Index", fitted="Fitted value")
  }
  if(missing(ylab)) ylab <- "Standardized residual"
  # ensure that horizontal grid line is drawn at 0
  breaks <- union(pretty(data[, "residual"]), 0)
  # create plot
  p <- ggplot(data) +
    geom_hline(aes(yintercept=-2.5), alpha=0.4) +
    geom_hline(aes(yintercept=2.5), alpha=0.4) +
    geom_point(mapping, size=size[1], ...)
  if(!is.null(labelData)) {
    # add labels for observations with largest distances
    labelMapping <- aes_string(x=abscissa, y="residual", label="index")
    p <- p + geom_text(labelMapping, data=labelData,
                       hjust=0, size=size[2], alpha=0.4)
  }
  p <- p + scale_y_continuous(breaks=breaks) +
    labs(title=main, x=xlab, y=ylab)
  if(!is.null(facets)) {
    # split plot into different panels
    if(length(facets) == 2) p <- p + facet_wrap(facets)
    else p <- p + facet_grid(facets)
  }
  p
}

# ----------------------

## plot robust distances (regression diagnostic plot)

rdiagPlot <- function(data, facets = attr(data, "facets"), size = c(2, 4),
                      id.n = NULL, main, xlab, ylab, ..., mapping) {
  ## initializations
  # extract data frame for vertical reference line
  lineData <- attr(data, "q")
  # check if robust distances are available
  msg <- "robust distances not available"
  by <- intersect(c("step", "fit"), names(data))
  if(length(by) > 0) {
    indices <- split(seq_len(nrow(data)), data[, by, drop=FALSE])
    onlyNA <- sapply(indices, function(i) all(is.na(data[i, "rd"])))
    if(all(onlyNA)) stop(msg)
    if(any(onlyNA)) {
      indices <- do.call(c, unname(indices[onlyNA]))
      data <- data[-indices, , drop=FALSE]
      lineData <- lineData[!onlyNA, , drop=FALSE]
      warning(msg, " for some submodels")
    }
  } else {
    onlyNA <- all(is.na(data[, "rd"]))
    if(onlyNA) stop(msg)
  }
  # define aesthetic mapping for regression diagnostic plot
  mapping <- aes_string(x="rd", y="residual", color="classification")
  ## construct data frame for labels
  labelData <- labelify(data, which="xyd", id.n=id.n)
  # define default title and axis labels
  if(missing(main)) main <- "Regression diagnostic plot"
  if(missing(xlab)) xlab <- "Robust distance computed by MCD"
  if(missing(ylab)) ylab <- "Standardized residual"
  # create plot
  p <- ggplot(data) +
    geom_hline(aes(yintercept=-2.5), alpha=0.4) +
    geom_hline(aes(yintercept=2.5), alpha=0.4)
  if(!is.null(lineData)) {
    # add reference line
    p <- p + geom_vline(aes_string(xintercept="q"), lineData, alpha=0.4)
  }
  p <- p + geom_point(mapping, size=size[1], ...)
  if(!is.null(labelData)) {
    # add labels for observations with largest distances
    labelMapping <- aes_string(x="rd", y="residual", label="index")
    p <- p + geom_text(labelMapping, data=labelData,
                       hjust=0, size=size[2], alpha=0.4)
  }
  p <- p + labs(title=main, x=xlab, y=ylab)
  if(!is.null(facets)) {
    # split plot into different panels
    if(length(facets) == 2) p <- p + facet_wrap(facets)
    else p <- p + facet_grid(facets)
  }
  p
}
