#' @include generics.R
{}



#' Default parallel coordinates plot.
#'
#' Code copied from function \code{\link[MASS]{parcoord}} of package
#' \code{MASS} to simply play arround with the visualization of archetypes.
#' At a later date, when it is clear which visualisation is the best, the
#' functionality is probabibly merged with the original function or it is
#' possible with parallel coordinate plots which are available et all.
#'
#' @param x A \eqn{n \times m} matrix or data frame who columns represent
#'   variables. Missing values are allowed.
#' @param col Line color.
#' @param lty Line type.
#' @param var.label Axes labels.
#' @param rx A \eqn{2 \times m} matrix with ranges for each dimension.
#' @param ... Passed to the underlying \code{\link[graphics]{matplot}}
#'   function.
#' @return Undefined.
#' @method pcplot default
#' @S3method pcplot default
pcplot.default <- function (x, col=gray(0.7), lty=1, var.label=TRUE,
                            rx=NULL, ...) {

  calcrx <- TRUE
  if ( is.null(rx) ) {
    rx <- apply(x, 2, range, na.rm=TRUE)
  }
  else {
    x <- rbind(rx, x)
    calcrx <- FALSE
  }

  sx <- sapply(1:ncol(x),
               function(i) {
                 (x[,i] - rx[1,i]) / (rx[2,i] - rx[1,i])
               })
  colnames(sx) <- colnames(x)

  x <- sx

  matplot(1:ncol(x), t(x), type="l", col=col, lty=lty,
          xlab="", ylab="", axes=FALSE, ...)
  axis(1, at=1:ncol(x), labels=colnames(x), ...)

  for (i in 1:ncol(x)) {
    lines(c(i, i), c(0, 1), col="grey70")
    if (var.label)
      text(c(i, i), c(0, 1), labels=format(rx[, i], digits=3),
           xpd=NA, offset=0.3, pos=c(1, 3), cex=0.7)
  }
}



#' Add lines to an existing parallel coordinates plot.
#' @param x A matrix or data frame containing the additional data.
#' @param data The data of the existing parallel coordinates plot.
#' @param col Line colors.
#' @param lty Line types.
#' @param ... Passed to underlying \code{\link[graphics]{matlines}}.
#' @return Undefined.
#' @rdname pcplot
lines.pcplot <- function(x, data, col=1, lty=1, ...) {
  rx <- apply(data, 2, range, na.rm=TRUE)

  x <- sapply(1:ncol(x),
              function(i) {
                (x[,i] - rx[1,i]) / (rx[2,i] - rx[1,i])
              })

  matlines(1:ncol(x), t(x), col=col, lty=lty, ...)
}
