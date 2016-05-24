#' smart boxplot
#'
#' draws points or boxes depending on sample size
#'
#' @param x a named list of numeric vectors
#' @param ... further arguments passed on to \code{\link{boxplot}}
#' @examples
#'
#' sboxplot(list(a=rnorm(15),b=rnorm(75),c=rnorm(10000)))
#'
sboxplot <- function(x, ...) {
  cuts <- c(20, 100)
  lens <- sapply(x, length)
  if (all(lens >= cuts[2])) {
    boxplot(x, range=0, ...)
  } else if (all(lens >= cuts[1] & lens < cuts[2])) {
    boxplot(x, ...)
  } else if (all(lens < cuts[1])) {
    stripchart(x, method="jitter", vertical=TRUE, pch=1, ...)
  } else {
    r <- range(unlist(x))
    d <- diff(r)
    n <- length(x)
    plot(0,0,type="n",xlab="",ylab="",xlim=c(0.5,n+.5),ylim=c(r[1] - .05 * d, r[2] + .05 * d),xaxt="n")
    for (i in seq_len(n)) {
      if (lens[i] < cuts[1]) {
        points(rep(i, lens[i]) + rnorm(lens[i], 0, .05), x[[i]])
      } else if (lens[i] >= cuts[1]) {
        box <- quantile(x[[i]], c(.25, .75))
        med <- median(x[[i]])
        polygon(i + c(-.1,-.1,.1,.1), c(box, rev(box)))
        segments(i - .1, med, i + .1, med, lwd=3)
        iqr <- diff(box)
        out <- box[1] - x[[i]] > 1.5 * iqr | x[[i]] - box[2] > 1.5 * iqr
        segments(i, box[1], i, min(x[[i]][!out]))
        segments(i, box[2], i, max(x[[i]][!out]))
        if (sum(out) > 0) {
          if (lens[i] < cuts[2]) {
            points(rep(i, sum(out)) + rnorm(sum(out), 0, .05), x[[i]][out])
          } else {
            segments(i, min(x[[i]][!out]), i, min(x[[i]]), lty=2)
            segments(i, max(x[[i]][!out]), i, max(x[[i]]), lty=2)
          }
        }
      }
    }
    axis(1, at=seq_len(n), labels=names(x))
  }
}
