window.bayescomm <-
function (x, start = NULL, end = NULL, thin = 1, ...) {
  oldstart <- x$call$start
  oldend <- x$call$its
  oldthin <- x$call$thin
  if (is.null(end)) {
    end <- oldend
  }
  if(is.null(start)) {
    start <- oldstart
  }
  if(is.null(thin)) {
    thin <- oldthin
  }
  trim <- function(x, start, end, thin, oldstart, oldthin) {
    suppressWarnings(as.matrix(window(x = mcmc(x,
                                               start = oldstart,
                                               thin = oldthin),
                                      start = start,
                                      end  = end,
                                      thin = thin)))
  }
  x$trace$R <- trim(x$trace$R, start, end, thin, oldstart, oldthin) 
  x$trace$z <- apply(x$trace$z, c(2, 3), trim, start, end, thin, oldstart, oldthin)
  x$trace$B <- lapply(x$trace$B, trim, start, end, thin, oldstart, oldthin)
  x$call$start <- start
  x$call$thin <- thin
  x$call$its <- end
  x
}

