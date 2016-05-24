approximation <- function(y, approx, n, type=c("standardized", "mean", "sum"),
                          approx.type=c("Edgeworth", "Saddlepoint")) {
  if (length(y) != length(approx)) {
    stop("Lengths of 'y' and 'approx' differ!")
  }
  type <- match.arg(type)
  approx.type <- match.arg(approx.type)
  value <- list(y=y, approx=approx, type=type, n=n, approx.type=approx.type)
  class(value) <- "approximation"
  return(value)
}

plot.approximation <- function(x, do.annotate=TRUE, ...) {
  extra.args <- list(...)
  default.args <- list(main=switch(x$type,
                         "standardized" = {
                           substitute(expression(paste(approx, "-Approximation of ",
                               Z==frac(S[n]-n%.%mu, sqrt(n*sigma^2)))),
                                      list(approx=x$approx.type))
                         },
                         "mean" = {
                           substitute(expression(paste(approx, "-Approximation of ",
                               Z==frac(1,n)%.%S[n])), list(approx=x$approx.type))
                         },
                         "sum" = {
                           substitute(expression(paste(approx, "-Approximation of ",
                               Z==S[n])), list(approx=x$approx.type)) 
                         }),                 
                       sub=expression(S[n]==Y[1]+cdots+Y[n]),
                       type="l",
                       xlab=expression(z),
                       ylab="Density")
  args <- c(list(x$y, x$approx), merge(extra.args, default.args))
  do.call("plot", args)
  if (do.annotate) {
    mtext(paste("n =", x$n), side=2, line=2)
  }
}

