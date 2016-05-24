#' Plot cv.enpls Object
#'
#' This function plots cv.enpls object.
#'
#' This function plots cv.enpls object.
#'
#' @param x An object of class \code{cv.enpls}.
#' @param main plot title
#' @param ... Other graphical parameters to be passed on to \code{plot}.
#'
#' @author Nan Xiao <\email{road2stat@@gmail.com}>
#'
#' @seealso See \code{\link{cv.enpls}} for performing ensemble PLS regression.
#'
#' @method plot cv.enpls
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
#' plot(cv.enpls.fit)

plot.cv.enpls = function(x, main = NULL, ...) {

  if (!inherits(x, 'cv.enpls'))
    stop('This function only works for objects of class "cv.enpls"')

  y.data = x$ypred
  y.real = y.data[, 'y.real']
  y.pred = y.data[, 'y.pred']

  plot(y.real, y.pred, 
       xlim = range(y.real), ylim = range(y.real), 
       xlab = 'Real Response', ylab = 'Predicted Response', 
       main = main, ...)
  abline(a = 0L, b = 1L)

}

#' Plot enpls.fs Object
#'
#' This function plots enpls.fs object.
#'
#' This function plots enpls.fs object.
#'
#' @param x An object of class \code{enpls.fs}.
#' @param sort Should the variables be sorted in decreasing order of importance?
#' @param nvar How many variables to show? Ignored if \code{sort = FALSE}.
#' @param main plot title
#' @param ... Other graphical parameters to be passed on to \code{dotchart}.
#'
#' @author Nan Xiao <\email{road2stat@@gmail.com}>
#'
#' @seealso See \code{\link{enpls.fs}} for feature selection with ensemble PLS.
#'
#' @method plot enpls.fs
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
#' plot(varimp)
#' plot(varimp, nvar = 10L)

plot.enpls.fs = function(x, 
                         sort = TRUE, nvar = NULL, 
                         main = NULL, ...) {

  if (!inherits(x, 'enpls.fs'))
    stop('This function only works for objects of class "enpls.fs"')

  varimp = x$variable.importance
  if (is.null(nvar)) nvar = length(varimp)

  if (sort == TRUE) {
    dotchart(sort(varimp, TRUE)[nvar:1], main = main, ...)
  } else {
    dotchart(rev(varimp), main = main, ...)
  }

}

#' Plot enpls.od Object
#'
#' This function plots enpls.od object.
#'
#' This function plots enpls.od object.
#'
#' @param x An object of class \code{enpls.od}.
#' @param criterion Criterion of being outlier, 
#' could be \code{'quantile'} or \code{'sd'}.
#' @param prob the quantile
#' @param sdtimes the times of sd
#' @param main plot title
#' @param ... Other graphical parameters to be passed on to \code{plot}.
#'
#' @author Nan Xiao <\email{road2stat@@gmail.com}>
#'
#' @seealso See \code{\link{enpls.od}} for outlier detection with ensemble PLS.
#'
#' @method plot enpls.od
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
#' plot(od, criterion = 'quantile')
#' plot(od, criterion = 'sd')

plot.enpls.od = function(x, 
                         criterion = c('quantile', 'sd'), 
                         prob = 0.05, sdtimes = 3L, 
                         main = NULL, ...) {

  if (!inherits(x, 'enpls.od'))
    stop('This function only works for objects of class "enpls.od"')

  criterion = match.arg(criterion)

  error.mean = x$error.mean
  error.sd = x$error.sd

  if (criterion == 'quantile') {
    vpos = quantile(error.mean, 1 - prob)
    hpos = quantile(error.sd, 1 - prob)
  } else {
    vpos = mean(error.mean) + (sdtimes * sd(error.mean))
    hpos = mean(error.sd) + (sdtimes * sd(error.sd))
  }

  yout = intersect(which(error.mean >= vpos), which(error.sd <= hpos))
  Xout = intersect(which(error.mean <= vpos), which(error.sd >= hpos))
  abnormal = intersect(which(error.mean >= vpos), which(error.sd >= hpos))

  plot(error.mean, error.sd, 
       xlab = 'Error Mean', ylab = 'Error SD', main = main, ...)

  abline(h = hpos, col = 'gray', lty = 2)
  abline(v = vpos, col = 'gray', lty = 2)

  if (length(yout) != 0L) text(error.mean[yout], error.sd[yout], 
                               labels = as.character(yout), 
                               col = 'red', cex = 0.7, pos = 3)
  if (length(Xout) != 0L) text(error.mean[Xout], error.sd[Xout], 
                               labels = as.character(Xout), 
                               col = 'blue', cex = 0.7, pos = 1)
  if (length(abnormal) != 0L) text(error.mean[abnormal], error.sd[abnormal], 
                                   labels = as.character(abnormal), 
                                   col = 'purple', cex = 0.7, pos = 1)

}
