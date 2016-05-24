plot.oemfit <- function(x, xvar = c("norm", "lambda", "loglambda", "dev"),
                        xlab = iname, ylab = "Coefficients",
                        ...) {
  xvar <- match.arg(xvar)
  nbeta <- as.matrix(x$beta)
  switch(xvar,
         "norm" = {
           index <- apply(abs(nbeta), 2, sum)
           iname <- "L1 Norm"
           xlim <- range(index)
         },
         "lambda" = {
           index <- x$lambda
           iname <- "Lambda"
           xlim <- rev(range(index))
         },
         "loglambda" = {
           index <- log(x$lambda)
           iname <- "Log Lambda"
           xlim <- rev(range(index))
         },
         "dev" = {
           index = x$sumSquare
           iname = "Sum of Squares"
           xlim <- range(index)
         }
       )
  matplot(index, t(nbeta), lty = 1, xlab = xlab, ylab = ylab, xlim = xlim,
          type = 'l', ...)
}
