#' Print cv.enpls Object
#'
#' This function prints cv.enpls object.
#'
#' This function prints cv.enpls object.
#'
#' @param x An object of class \code{cv.enpls}.
#' @param ... Other parameters to be passed on to \code{print}.
#'
#' @author Nan Xiao <\email{road2stat@@gmail.com}>
#'
#' @seealso See \code{\link{cv.enpls}} for performing ensemble PLS regression.
#'
#' @method print cv.enpls
#'
#' @export
#'
#' @examples
#' data(alkanes)
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' cv.enpls.fit = cv.enpls(x, y, MCtimes = 20)
#' print(cv.enpls.fit)

print.cv.enpls = function(x, ...) {

  if (!inherits(x, 'cv.enpls'))
    stop('This function only works for objects of class "cv.enpls"')

  cat('Cross Validation Result for Ensemble Partial Least Squares\n')
  cat('---\n')
  cat(paste('RMSE = ', sprintf("%.4f", x$RMSE), 
            ', R2 = ', sprintf("%.6f", x$R2), '\n', 
            sep = ''))

}

#' Print Fitted Ensemble Partial Least Squares Object
#'
#' This function prints coefficients of each model in the enpls.en object.
#'
#' This function prints coefficients of each model in the enpls.en object.
#'
#' @param x An object of class \code{enpls.en}.
#' @param ... Other parameters to be passed on to \code{print}.
#'
#' @author Nan Xiao <\email{road2stat@@gmail.com}>
#'
#' @seealso See \code{\link{enpls.en}} for performing ensemble PLS regression.
#'
#' @method print enpls.en
#'
#' @export
#'
#' @examples
#' data(alkanes)
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' enpls.fit = enpls.en(x, y, MCtimes = 100)
#' print(enpls.fit)

print.enpls.en = function(x, ...) {

  if (!inherits(x, 'enpls.en'))
    stop('This function only works for objects of class "enpls.en"')

  coefmeta = coef(x[[1]][[1]], intercept = TRUE)[, 1, 1]
  varcount = length(coefmeta)
  mctimes  = length(x)
  coefdf   = matrix(NA, ncol = mctimes, nrow = varcount)
  for (i in 1:mctimes) coefdf[, i] = coef(x[[i]][[1]], intercept = TRUE)[, 1, 1]
  rownames(coefdf) = names(coefmeta)

  cat('Coefficients of the Models by Ensemble Partial Least Squares\n')
  cat('---\n')
  print(coefdf)

}

#' Print enpls.fs Object
#'
#' This function prints enpls.fs object.
#'
#' This function prints enpls.fs object.
#'
#' @param x An object of class \code{enpls.fs}.
#' @param sort Should the variables be sorted in decreasing order of importance?
#' @param nvar How many variables to show? Ignored if \code{sort = FALSE}.
#' @param ... Other parameters to be passed on to \code{print}.
#'
#' @author Nan Xiao <\email{road2stat@@gmail.com}>
#'
#' @seealso See \code{\link{enpls.fs}} for feature selection with ensemble PLS.
#'
#' @method print enpls.fs
#'
#' @export
#'
#' @examples
#' data(alkanes)
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' varimp = enpls.fs(x, y, MCtimes = 100)
#' print(varimp)
#' print(varimp, nvar = 10L)

print.enpls.fs = function(x, sort = TRUE, nvar = NULL, ...) {

  if (!inherits(x, 'enpls.fs'))
    stop('This function only works for objects of class "enpls.fs"')

  varimp = x$variable.importance
  if (is.null(nvar)) nvar = length(varimp)

  cat('Variable Importance by Ensemble Partial Least Squares\n')
  cat('---\n')
  if (sort == TRUE) {
    print(data.frame('Importance' = sort(varimp, TRUE)[1:nvar]))
  } else {
    print(data.frame('Importance' = varimp))
  }

}

#' Print enpls.od Object
#'
#' This function prints enpls.od object.
#'
#' This function prints enpls.od object.
#'
#' @param x An object of class \code{enpls.od}.
#' @param ... Other parameters to be passed on to \code{print}.
#'
#' @author Nan Xiao <\email{road2stat@@gmail.com}>
#'
#' @seealso See \code{\link{enpls.od}} for outlier detection with ensemble PLS.
#'
#' @method print enpls.od
#'
#' @export
#'
#' @examples
#' data(alkanes)
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' od = enpls.od(x, y, MCtimes = 100)
#' print(od)

print.enpls.od = function(x, ...) {

  if (!inherits(x, 'enpls.od'))
    stop('This function only works for objects of class "enpls.od"')

  cat('Outlier Detection by Ensemble Partial Least Squares\n')
  cat('---\n')
  cat('Mean residual for each sample:\n')
  print(x$error.mean)
  cat('---\n')
  cat('Residual SD for each sample:\n')
  print(x$error.sd)

}
