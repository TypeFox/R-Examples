#' An S4 class that stores the outputs of the fitted IGLS model.
#'
#' An MLwiN model run via the IGLS estimation method is represented by an "mlwinfitIGLS" object
#'
#' @section An instance of the Class:
#'  An instance is created by calling function \code{\link{runMLwiN}}.
#'
#' @slot Nobs Computes the number of complete observations.
#' @slot DataLength Total number of cases.
#' @slot Hierarchy For each higher level of a multilevel model, returns the number of units at that level, together with the minimum, mean and maximum number of lower-level units nested within units of the current level.
#' @slot D A vector specifying the type of distribution to be modelled, which can include \code{'Normal'}, \code{'Binomial'} \code{'Poisson'}, \code{'Multinomial'}, \code{'Multivariate Normal'}, or \code{'Mixed'}.
#' @slot Formula A formula object (or a character string) specifying a multilevel model.
#' @slot levID A character string (vector) of the specified level ID(s).
#' @slot FP Displays the fixed part estimates.
#' @slot RP Displays the random part estimates.
#' @slot FP.cov Displays a covariance matrix of the fixed part estimates.
#' @slot RP.cov Displays a covariance matrix of the random part estimates.
#' @slot elapsed.time Calculates the CPU time used for fitting the model.
#' @slot call The matched call.
#' @slot LIKE The deviance statistic (-2*log(like)).
#' @slot Converged Boolean indicating whether the model has converged
#' @slot Iterations Number of iterations that the model has run for
#' @slot Meth If \code{Meth = 0} estimation method is set to RIGLS. If \code{Meth = 1} estimation method is set to IGLS.
#' @slot residual If \code{resi.store} is \code{TRUE}, then the residual estimates at all levels are returned.
#' @slot data The data.frame that was used to fit the model.
#' @slot nonlinear A character vector specifying linearisation method used. The first element specifies marginal quasi-likelihood linearization (\code{N = 0}) or penalised quasi-likelihood linearization (\code{N = 1}); The second element specifies first (\code{M = 1}) or second (\code{M = 2}) order approximation.
#' @slot version The MLwiN version used to fit the model
#'
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#'
#' @seealso
#' \code{\link{runMLwiN}}
#'
#' @examples
#' \dontrun{
#' library(R2MLwiN)
#' # NOTE: if MLwiN not saved in location R2MLwiN defaults to, specify path via:
#' # options(MLwiN_path = 'path/to/MLwiN vX.XX/')
#' # If using R2MLwiN via WINE, the path may look like this:
#' # options(MLwiN_path = '/home/USERNAME/.wine/drive_c/Program Files (x86)/MLwiN vX.XX/')
#'
#' ## Example: tutorial
#' data(tutorial, package = "R2MLwiN")
#'
#' (mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student),
#'                      data = tutorial))
#'
#' ##summary method
#' summary(mymodel)
#'
#' ##logLik method
#' logLik(mymodel)
#' }
#'
#' @name mlwinfitIGLS-class
#' @rdname mlwinfitIGLS-class
#' @exportClass mlwinfitIGLS
setClass(Class = "mlwinfitIGLS", representation = representation(version = "character", Nobs = "numeric", DataLength = "numeric",
                                                                 Hierarchy = "ANY", D = "ANY", Formula = "ANY", levID = "character", FP = "numeric", RP = "numeric", RP.cov = "matrix",
                                                                 FP.cov = "matrix", LIKE = "ANY", elapsed.time = "numeric", call = "ANY", residual = "list", Converged = "logical",
                                                                 Iterations = "numeric", Meth = "numeric", nonlinear = "numeric", data = "data.frame"))

#' Extract or Replace parts of "mlwinfitIGLS" objects
#' @param x data frame
#' @param i,j elements to extract or replace. For \code{[} and \code{[[}, these are \code{character}.
#' @param drop not used.
#' @param value a suitable replacement value.
#' @rdname extract-methods-igls
setMethod("[", "mlwinfitIGLS", function(x, i, j, drop) {
  if (i == "version") {
    return(x@version)
  }
  if (i == "Nobs") {
    return(x@Nobs)
  }
  if (i == "DataLength") {
    return(x@DataLength)
  }
  if (i == "Hierarchy") {
    return(x@Hierarchy)
  }
  if (i == "D") {
    return(x@D)
  }
  if (i == "Formula") {
    return(x@Formula)
  }
  if (i == "levID") {
    return(x@levID)
  }
  if (i == "FP") {
    return(x@FP)
  }
  if (i == "RP") {
    return(x@RP)
  }
  if (i == "FP.cov") {
    return(x@FP.cov)
  }
  if (i == "RP.cov") {
    return(x@RP.cov)
  }
  if (i == "elapsed.time") {
    return(x@elapsed.time)
  }
  if (i == "call") {
    return(x@call)
  }
  if (i == "LIKE") {
    return(x@LIKE)
  }
  if (i == "Converged") {
    return(x@Converged)
  }
  if (i == "Iterations") {
    return(x@Iterations)
  }
  if (i == "Meth") {
    return(x@Meth)
  }
  if (i == "nonlinear") {
    return(x@nonlinear)
  }
  if (i == "residual") {
    return(x@residual)
  }
  if (i == "data") {
    return(x@data)
  }
})

#' @rdname extract-methods-igls
setReplaceMethod("[", signature(x = "mlwinfitIGLS"), function(x, i, j, value) {
  if (i == "version") {
    x@version <- value
  }
  if (i == "Nobs") {
    x@Nobs <- value
  }
  if (i == "DataLength") {
    x@DataLength <- value
  }
  if (i == "Hierarchy") {
    x@Hierarchy <- value
  }
  if (i == "D") {
    x@D <- value
  }
  if (i == "Formula") {
    x@Formula <- value
  }
  if (i == "levID") {
    x@levID <- value
  }
  if (i == "FP") {
    x@FP <- value
  }
  if (i == "RP") {
    x@RP <- value
  }
  if (i == "FP.cov") {
    x@FP.cov <- value
  }
  if (i == "RP.cov") {
    x@RP.cov <- value
  }
  if (i == "elapsed.time") {
    x@elapsed.time <- value
  }
  if (i == "call") {
    x@call <- value
  }
  if (i == "LIKE") {
    x@LIKE <- value
  }
  if (i == "Converged") {
    x@Converged <- value
  }
  if (i == "Iterations") {
    x@Iterations <- value
  }
  if (i == "Meth") {
    x@Meth <- value
  }
  if (i == "nonlinear") {
    x@nonlinear <- value
  }
  if (i == "residual") {
    x@residual <- value
  }
  if (i == "data") {
    x@data <- value
  }
  validObject(x)
  return(x)
})

#' @rdname extract-methods-igls
setMethod("[[", "mlwinfitIGLS", function(x, i, j, drop) {
  if (i == "version") {
    return(x@version)
  }
  if (i == "Nobs") {
    return(x@Nobs)
  }
  if (i == "DataLength") {
    return(x@DataLength)
  }
  if (i == "Hierarchy") {
    return(x@Hierarchy)
  }
  if (i == "D") {
    return(x@D)
  }
  if (i == "Formula") {
    return(x@Formula)
  }
  if (i == "levID") {
    return(x@levID)
  }
  if (i == "FP") {
    return(x@FP)
  }
  if (i == "RP") {
    return(x@RP)
  }
  if (i == "FP.cov") {
    return(x@FP.cov)
  }
  if (i == "RP.cov") {
    return(x@RP.cov)
  }
  if (i == "elapsed.time") {
    return(x@elapsed.time)
  }
  if (i == "call") {
    return(x@call)
  }
  if (i == "LIKE") {
    return(x@LIKE)
  }
  if (i == "Converged") {
    return(x@Converged)
  }
  if (i == "Iterations") {
    return(x@Iterations)
  }
  if (i == "Meth") {
    return(x@Meth)
  }
  if (i == "nonlinear") {
    return(x@nonlinear)
  }
  if (i == "residual") {
    return(x@residual)
  }
  if (i == "data") {
    return(x@data)
  }
})

#' @rdname extract-methods-igls
setReplaceMethod("[[", signature(x = "mlwinfitIGLS"), function(x, i, j, value) {
  if (i == "version") {
    x@version <- value
  }
  if (i == "Nobs") {
    x@Nobs <- value
  }
  if (i == "DataLength") {
    x@DataLength <- value
  }
  if (i == "Hierarchy") {
    x@Hierarchy <- value
  }
  if (i == "D") {
    x@D <- value
  }
  if (i == "Formula") {
    x@Formula <- value
  }
  if (i == "levID") {
    x@levID <- value
  }
  if (i == "FP") {
    x@FP <- value
  }
  if (i == "RP") {
    x@RP <- value
  }
  if (i == "FP.cov") {
    x@FP.cov <- value
  }
  if (i == "RP.cov") {
    x@RP.cov <- value
  }
  if (i == "elapsed.time") {
    x@elapsed.time <- value
  }
  if (i == "call") {
    x@call <- value
  }
  if (i == "LIKE") {
    x@LIKE <- value
  }
  if (i == "Converged") {
    x@Converged <- value
  }
  if (i == "Iterations") {
    x@Iterations <- value
  }
  if (i == "Meth") {
    x@Meth <- value
  }
  if (i == "nonlinear") {
    x@nonlinear <- value
  }
  if (i == "residual") {
    x@residual <- value
  }
  if (i == "data") {
    x@data <- value
  }
  validObject(x)
  return(x)
})

#' Summarize "mlwinfitIGLS" objects
#' @param object,x an \code{\link{mlwinfitIGLS-class}} object
#' @param ... other parameters
#' @param digits the number of significant digits to use when printing.
#' @param signif.stars logical. If TRUE, 'significance stars' are printed for each coefficient.
#' @rdname summary-methods-igls
#' @export
setMethod("summary", signature(object = "mlwinfitIGLS"), function(object, ...) {
  object
})

printIGLS <- function(x, digits = max(3, getOption("digits") - 2), signif.stars = getOption("show.signif.stars"),
                      ...) {
  object <- summary(x)
  align2right <- function(titlename, ele) {
    # for printing the table on the screen
    all.ele <- c(titlename, ele)
    len.all.ele <- nchar(all.ele)
    max.len.ele <- max(len.all.ele)
    for (j in 1:length(all.ele)) {
      if (len.all.ele[j] < max.len.ele) {
        len.diff <- max.len.ele - len.all.ele[j]
        all.ele[j] <- paste(paste(rep(" ", len.diff), collapse = ""), all.ele[j], sep = "")
      }
    }
    all.ele
  }
  
  align2left <- function(titlename, ele) {
    # for printing the table on the screen
    all.ele <- c(titlename, ele)
    len.all.ele <- nchar(all.ele)
    max.len.ele <- max(len.all.ele)
    for (j in 1:length(all.ele)) {
      if (len.all.ele[j] < max.len.ele) {
        len.diff <- max.len.ele - len.all.ele[j]
        all.ele[j] <- paste(all.ele[j], paste(rep(" ", len.diff), collapse = ""), sep = "")
      }
    }
    all.ele
  }
  
  signifstar <- function(pval) {
    starstr <- "N/A"
    if (!is.na(pval) && pval >= 0 && pval <= 1) {
      if (pval < 0.001) {
        starstr <- "***"
      }
      if (pval >= 0.001 && pval < 0.01) {
        starstr <- "** "
      }
      if (pval >= 0.01 && pval < 0.05) {
        starstr <- "*  "
      }
      if (pval >= 0.05 && pval < 0.1) {
        starstr <- ".  "
      }
      if (pval >= 0.1) {
        starstr <- "   "
      }
    }
    starstr
  }
  cat("\n")
  cat(paste(rep("-", 50), collapse = "*"), "\n")
  cat(object@version, " multilevel model", paste("(", object@D[1], ")", sep = ""), "\n")
  if (!is.null(object@Hierarchy))
    print(object@Hierarchy)
  cat("Estimation algorithm:  ")
  if (object@Meth == 1) {
    cat("IGLS")
  } else {
    cat("RIGLS")
  }
  if (object@D[1] != "Normal" && object@D[1] != "Multivariate Normal") {
    if (object@nonlinear[1] == 0) {
      cat(" MQL")
    }
    if (object@nonlinear[1] == 1) {
      cat(" PQL")
    }
    cat(object@nonlinear[2])
  }
  cat("        Elapsed time :", paste(round(object@elapsed.time, 2), "s", sep = ""), "\n")
  cat("Number of obs: ", object@Nobs, paste0("(from total ", object@DataLength, ")"))
  if (object@Converged) {
    cat("        The model converged after", object@Iterations, "iterations.\n")
  } else {
    cat("        The model did not converge after", object@Iterations, "iterations.\n")
  }
  cat(paste("Log likelihood:     ", round(-0.5 * object@LIKE, 1)), "\n")
  cat(paste("Deviance statistic: ", round(object@LIKE, 1)), "\n")
  cat(paste(rep("-", 50), collapse = "-"), "\n")
  cat("The model formula:\n")
  print(object@Formula)  #cat(gsub('[[:space:]]','',object@Formula),'\n')
  levID0 <- object@levID
  levID.display <- ""
  if (is.na(levID0[length(levID0)])) {
    levID0 <- levID0[-length(levID0)]
  }
  for (i in 1:length(levID0)) {
    levID.display <- paste(levID.display, "Level ", length(levID0) + 1 - i, ": ", levID0[i], "     ", sep = "")
  }
  cat(levID.display, "\n")
  cat(paste(rep("-", 50), collapse = "-"), "\n")
  
  FP.names <- names(object@FP)
  RP.names <- names(object@RP)
  
  cat("The fixed part estimates: ", "\n")
  FP.print <- rbind(object@FP, sqrt(diag(object@FP.cov)))
  z.score <- FP.print[1, ]/FP.print[2, ]
  p.value <- 2 * pnorm(abs(z.score), lower.tail = FALSE)
  strstar <- as.vector(sapply(p.value, signifstar))
  qt025 <- FP.print[1, ] - qnorm(0.975) * FP.print[2, ]
  qt975 <- FP.print[1, ] + qnorm(0.975) * FP.print[2, ]
  FP.print <- rbind(FP.print, z.score, p.value, qt025, qt975)
  FP.names2 <- gsub("FP+\\_", "", FP.names)
  
  printcol0 <- align2left("        ", FP.names2)
  printcol1 <- align2right("Coef.", format(round(FP.print[1, ], digits), nsmall = digits))
  printcol2 <- align2right("Std. Err.", format(round(FP.print[2, ], digits), nsmall = digits))
  printcol3 <- align2right("z", format(round(FP.print[3, ], 2), nsmall = 2))
  printcol4 <- align2right("Pr(>|z|)", formatC(FP.print[4, ]))
  printcol4b <- align2right("   ", strstar)
  printcol5 <- align2right("[95% Conf.", format(round(FP.print[5, ], digits), nsmall = digits))
  printcol6 <- align2right("Interval]", format(round(FP.print[6, ], digits), nsmall = digits))
  for (i in 1:(ncol(FP.print) + 1)) {
    cat(printcol0[i], " ", printcol1[i], " ", printcol2[i], " ", printcol3[i], " ", printcol4[i], " ", printcol4b[i],
        " ", printcol5[i], " ", printcol6[i], "\n")
  }
  if (signif.stars) {
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ", "\n")
  }
  nlev <- length(object@levID)
  if (is.na(object@levID[length(object@levID)])) {
    mlwinlev <- (nlev - 1):1
    levID2 <- levID0
  } else {
    mlwinlev <- nlev:1
    levID2 <- object@levID
  }
  
  RP.print <- rbind(object@RP, sqrt(diag(object@RP.cov)))
  qt025 <- RP.print[1, ] - qnorm(0.975) * RP.print[2, ]
  qt975 <- RP.print[1, ] + qnorm(0.975) * RP.print[2, ]
  RP.print <- rbind(RP.print, qt025, qt975)
  for (i in 1:length(mlwinlev)) {
    RPx.pos <- grep(paste("RP", mlwinlev[i], sep = ""), RP.names)
    if (length(RPx.pos) != 0) {
      cat(paste(rep("-", 50), collapse = "-"), "\n")
      RPx.names <- gsub(paste("RP+", mlwinlev[i], "+\\_", sep = ""), "", RP.names[RPx.pos])
      RPx <- as.matrix(RP.print[, RPx.pos], nrow = 4)
      printcol0 <- align2left("        ", RPx.names)
      printcol1 <- align2right("Coef.", format(round(RPx[1, ], digits), nsmall = digits))
      printcol2 <- align2right("Std. Err.", format(round(RPx[2, ], digits), nsmall = digits))
      cat("The random part estimates at the", levID2[i], "level:", "\n")
      for (i in 1:(ncol(RPx) + 1)) {
        cat(printcol0[i], " ", printcol1[i], " ", printcol2[i], "\n")
      }
    }
  }
  cat(paste(rep("-", 50), collapse = "*"), "\n")
}

#' @rdname summary-methods-igls
#' @export
setMethod("print", "mlwinfitIGLS", printIGLS)

#' @rdname summary-methods-igls
#' @export
setMethod("show", signature(object = "mlwinfitIGLS"), function(object) printIGLS(object))

updateMLwiN <- function(object, Formula., levID., estoptions., ..., keep.order = TRUE, evaluate = TRUE) {
  update.formula2 <- function (old, new, ...) 
  {
      C_updateform <- get("C_updateform", asNamespace("stats"), inherits=FALSE)
      tmp <- .Call(C_updateform, as.formula(old), as.formula(new))
      out <- formula(terms.formula(tmp, simplify = FALSE))
      return(out)
  }
  environment(update.formula2) <- environment(update.formula)
  my.update.formula <- function(old, new, keep.order = TRUE, ...) {
    env <- environment(as.formula(old))
    tmp <- update.formula2(as.formula(old), as.formula(new))
    out <- formula(terms.formula(tmp, simplify = FALSE, keep.order = keep.order))
    environment(out) <- env
    return(out)
  }
  if (is.null(newcall <- getCall(object)))
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (length(newcall$Formula))
    newcall$Formula <- eval(newcall$Formula)
  if (!missing(Formula.)) {
    newcall$Formula <- my.update.formula(as.formula(newcall$Formula), Formula., keep.order = keep.order)
  }
  if (!missing(levID.)) {
    newcall$levID <- {
      if (length(newcall$levID))
        my.update.formula(as.formula(newcall$levID), levID., keep.order = keep.order) else levID.
    }
  }
  if (!missing(estoptions.)) {
    newcall$estoptions <- {
      if (length(newcall$estoptions))
        my.update.formula(as.formula(newcall$estoptions), estoptions., keep.order = keep.order) else estoptions.
    }
  }
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(newcall)))
    for (a in names(extras)[existing]) newcall[[a]] <- extras[[a]]
    if (any(!existing)) {
      newcall <- c(as.list(newcall), extras[!existing])
      newcall <- as.call(newcall)
    }
  }
  if (evaluate)
    eval(newcall, sys.parent()) else newcall
}

#' Update "mlwinfitIGLS" objects
#' @param object a valid \code{mlwinfitIGLS} class object with an R function call component named \code{call}, the expression used to create itself.
#' @param Formula. changes to the formula. This is a two sided formula where "." is substituted for existing components in the \code{Formula} component of \code{object$call}.
#' @param levID. changes to the specifications of level ID(s).
#' @param estoptions. changes to the specifications of a list of options used for estimating the model.
#' @param ...  additional arguments to the call, or arguments with changed values.
#' @param keep.order a logical value indicating whether the terms should keep their positions.
#' @param evaluate  if \code{TRUE} (the default) the new call is evaluated;
#' otherwise the call is returned as an unevaluated expression.
#' @return either a new updated \code{mlwinfitIGLS} class object, or else an unevaluated expression for creating such an object.
#' @export
setMethod("update", signature(object = "mlwinfitIGLS"), updateMLwiN)

#' "mlwinfitIGLS" model formula
#' @param x See \code{\link[stats]{formula}}
#' @param env See \code{\link[stats]{formula}}
#' @param ... Other arguments; see \code{\link[stats]{formula}}
#' @export
setMethod("formula", "mlwinfitIGLS", function(x, env = parent.frame(), ...) {
  as.formula(x@Formula)
})

#' Extract the coefficient vector from "mlwinfitIGLS" objects
#' @param object An \code{\link{mlwinfitIGLS-class}} object
#' @param ... Other arguments
#' @seealso \code{\link[stats]{coef}}
#' @export
setMethod("coef", signature(object = "mlwinfitIGLS"), function(object, ...) {
  c(object@FP, object@RP)
})

#' @rdname coef-mlwinfitIGLS-method
#' @export
setMethod("coefficients", signature(object = "mlwinfitIGLS"), function(object, ...) {
  coef(object)
})

#' Extract the approximate variance-covariance matrix from "mlwinfitIGLS" objects
#' @param object An \code{\link{mlwinfitIGLS-class}} object
#' @param ... Other arguments
#' @seealso \code{\link[stats]{vcov}}
#' @export
setMethod("vcov", signature(object = "mlwinfitIGLS"), function(object, ...) {
  m <- matrix(0, nrow(object@FP.cov) + nrow(object@RP.cov), ncol(object@FP.cov) + ncol(object@RP.cov))
  colnames(m) <- c(colnames(object@FP.cov), colnames(object@RP.cov))
  rownames(m) <- c(rownames(object@FP.cov), rownames(object@RP.cov))
  m[colnames(object@FP.cov), rownames(object@FP.cov)] <- object@FP.cov
  m[colnames(object@RP.cov), rownames(object@RP.cov)] <- object@RP.cov
  m
})

#' Returns the residual degrees-of-freedom extracted from "mlwinfitIGLS" objects.
#' @param object An \code{\link{mlwinfitIGLS-class}} object.
#' @param ... Other arguments
#' @seealso \code{\link[stats]{nobs}}, \code{\link[stats]{coef}}
#' @export
setMethod("df.residual", signature(object = "mlwinfitIGLS"), function(object, ...) {
  nobs(object) - length(coef(object))
})

#' Returns the fitted values from "mlwinfitIGLS" objects.
#' @param object An \code{\link{mlwinfitIGLS-class}} object.
#' @param ... Other arguments.
#' @seealso \code{\link[stats]{fitted}}
#' @export
setMethod("fitted", signature(object = "mlwinfitIGLS"), function(object, ...) {
  predict(object, type = "response")
})

#' @rdname fitted-mlwinfitIGLS-method
#' @export
setMethod("fitted.values", signature(object = "mlwinfitIGLS"), function(object, ...) {
  fitted(object)
})

#' Returns the residual data from "mlwinfitIGLS" objects.
#'
#' @param object An \code{\link{mlwinfitIGLS-class}} object
#' @param ... Other arguments.
#' @seealso \code{\link[stats]{residuals}}
#' @export
setMethod("residuals", signature(object = "mlwinfitIGLS"), function(object, ...) {
  form <- Formula.translate(object@Formula, object@D, object@data)
  if (!is.list(form$resp)) {
    D <- object@D
    indata <- object@data
    tval <- fitted(object)
    if (D[1] == "Binomial") {
      tval <- tval * indata[, D[3]]
    }
    if (D[1] == "Poisson" || D[1] == "Negbinom") {
      if (!is.na(D[3])) {
        tval <- tval + indata[, D[3]]
      }
    }
    object@data[[form$resp]] - tval
  } else {
    warning("residuals only implemented for univariate models")
    NULL
  }
})

#' @rdname residuals-mlwinfitIGLS-method
#' @export
setMethod("resid", signature(object = "mlwinfitIGLS"), function(object, ...) {
  residuals(object)
})

#' Returns the predicted data from "mlwinfitIGLS" objects.
#' @param object An \code{\link{mlwinfitIGLS-class}} object.
#' @param newdata data frame for which to evaluate predictions
#' @param params a character vector specifying the parameters to use in evaluating predictions.
#' If \code{NULL}, \code{names(object[["FP"]])} is used by default.
#' @param type when this has the value \code{"link"} (default) the linear predictor is returned.
#' When \code{type="terms"} each component of the linear predictor is returned seperately. When \code{type="response"} predictions on the scale of the response are returned.
#' @param se.fit logical. When this is \code{TRUE} (not default) standard error estimates are returned for each prediction.
#' @param terms if \code{type="terms"}, which terms (default is all terms), a character vector.
#' @param ... Other arguments.
#' @seealso \code{\link[stats]{predict}}
#' @export
setMethod("predict", signature(object = "mlwinfitIGLS"), function(object, newdata = NULL, params = NULL, type = "link", se.fit = FALSE,
                                              terms = NULL, ...) {
  if (is.null(newdata)) {
    indata <- object@data
  } else {
    indata <- Formula.translate(object@Formula, object@D, object@data)$indata
  }
  if (is.null(params)) {
    fp.names <- names(FP <- object@FP)
  } else {
    fp.names <- params
  }
  x.names <- sub("FP_", "", fp.names)
  if (type == "link") {
    tval <- as.vector(as.matrix(indata[x.names]) %*% as.matrix(object@FP[fp.names])[, 1])
    if (se.fit) {
      # seval <- as.vector(sqrt(diag(as.matrix(indata[x.names]) %*% as.matrix(object@FP.cov[fp.names, fp.names]) %*%
      # t(as.matrix(indata[x.names])))))
      seval <- as.vector(sqrt(rowSums(as.matrix(indata[x.names]) %*% as.matrix(object@FP.cov[fp.names, fp.names]) *
                                        indata[x.names])))
      return(list(fit = tval, se.fit = seval))
    } else {
      return(tval)
    }
  }
  if (type == "terms") {
    if (!is.null(terms)) {
      x.names <- terms
      fp.names <- paste0("FP_", terms)
    }
    tval <- as.matrix(t(t(indata[x.names]) * object@FP[fp.names]))
    if (se.fit) {
      seval <- as.matrix(sqrt(t(t(indata[x.names]^2) * diag(object@FP.cov[fp.names, fp.names]))))
      return(list(fit = tval, se.fit = seval))
    } else {
      return(tval)
    }
  }
  if (type == "response") {
    tval <- as.vector(as.matrix(indata[x.names]) %*% as.matrix(object@FP[fp.names])[, 1])
    names(tval) <- 1:(length(tval))
    D <- object@D
    if (D[1] == "Normal" || D[1] == "Multivariate Normal") {
      return(tval)
    }
    if (D[1] == "Binomial") {
      if (D[2] == "logit") {
        antilogit <- function(x) {
          exp(x)/(1 + exp(x))
        }
        return(antilogit(tval) * indata[, D[3]])
      }
      if (D[2] == "probit") {
        return(pnorm(tval) * indata[, D[3]])
      }
      if (D[2] == "cloglog") {
        anticloglog <- function(x) {
          1 - exp(-exp(x))
        }
        return(anticloglog(tval) * indata[, D[3]])
      }
    }
    if (D[1] == "Poisson") {
      if (is.na(D[3])) {
        return(exp(tval))
      } else {
        return(exp(tval + indata[, D[3]]))
      }
    }
    if (D[1] == "Negbinom") {
      if (is.na(D[3])) {
        return(exp(tval))
      } else {
        return(exp(tval + indata[, D[3]]))
      }
    }
    if (D[1] == "Mixed") {
    }
    if (D[1] == "Multinomial") {
    }
    warning("link function transformation not yet implemented")
    return(NULL)
  }
})

#' Returns the log-likelihood from "mlwinfitIGLS" objects.
#' @param object An \code{\link{mlwinfitIGLS-class}} object.
#' @param ... Other arguments.
#' @seealso \code{\link[stats]{logLik}}
#' @export
setMethod("logLik", signature(object = "mlwinfitIGLS"), function(object, ...) {
  D <- object@D
  if (D[1] == "Normal" || D[1] == "Multivariate Normal") {
    val <- -0.5 * deviance(object)
    attr(val, "df") <- length(coef(object))
    attr(val, "nobs") <- nobs(object)
    class(val) <- "logLik"
    return(val)
  } else {
    warning("logLik only available for Normal models")
    return(NULL)
  }
})

#' Returns the deviance from "mlwinfitIGLS" objects.
#' @param object An \code{\link{mlwinfitIGLS-class}} object
#' @param ... Other arguments
#' @seealso \code{\link[stats]{deviance}}
#' @export
setMethod("deviance", signature(object = "mlwinfitIGLS"), function(object, ...) {
  D <- object@D
  if (D[1] == "Normal" || D[1] == "Multivariate Normal") {
    return(object@LIKE)
  } else {
    warning("deviance only available for Normal models")
    return(NULL)
  }
})

#' Returns the number of used observations from "mlwinfitIGLS" objects.
#' @param object An \code{\link{mlwinfitIGLS-class}} object.
#' @param ... Other arguments.
#' @seealso \code{\link[stats]{nobs}}
#' @export
setMethod("nobs", signature(object = "mlwinfitIGLS"), function(object, ...) {
  object@Nobs
})
