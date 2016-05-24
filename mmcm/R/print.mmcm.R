#' Print function for mmcm object
#' 
#' This function print result of function \code{\link{mcm.mvt}}, \code{\link{mmcm.mvt}} and \code{\link{mmcm.resamp}}
#' 
#' The case where printed "More than 2 contrast coefficient vectors were selected", some contrast may be unsuitable.
#' 
#' @param x Object of class \code{mmcm}, which is result of function \code{\link{mcm.mvt}}, \code{\link{mmcm.mvt}} and \code{\link{mmcm.resamp}}.
#' @param digits a non-null value for digits specifies the minimum number of significant digits to be printed in values. The default, NULL, uses \code{\link[base:options]{getOption}}(digits). (For the interpretation for complex numbers see \code{\link[base:Round]{signif}}.) Non-integer values will be rounded down, and only values greater than or equal to 1 and no greater than 22 are accepted.
#' @param ... Further arguments passed to or from other methods.
#' @method print mmcm
#' @seealso
#' \code{\link[base:print]{print.default}},
#' \code{\link{mmcm.mvt}},
#' \code{\link{mmcm.resamp}},
#' \code{\link{mcm.mvt}}
#' @keywords print
#' @export
print.mmcm <- function(x, digits = getOption("digits"), ...) {

  ####################
  # print result
  ####################
  if (x$contrast == "More than 2 contrast coefficient vectors were selected") {
    class(x) <- "htest"
    print(x)
    cat("More than 2 contrast coefficient vectors were selected\n\n")
  } else {
    class(x) <- "htest"
    print(x)
    msg <- paste(
      names(x$contrast      ), " = ", x$contrast      , ", ",
      names(x$contrast.index), " = ", x$contrast.index, "\n\n",
      names(x$error         ), " = ", format.pval(x$error, digits = digits), "\n",
      names(x$msg           ), " = ", x$msg           , "\n\n",
      sep = ""
    )
    cat(msg)
  }
  invisible(0)
  
}

