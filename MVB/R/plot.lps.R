plot.lps <- function(x, xvar = c("norm", "lambda", "loglambda", "dev"),
                     xlab = iname, ylab = "Coefficients",
                        ...) {
  xvar <- match.arg(xvar)
  nbeta <- as.matrix(x$beta)
  switch(xvar,
         "norm" = {
           index <- apply(abs(nbeta), 2, sum)
           iname <- "L1 Norm"
         },
         "lambda" = {
           index <- -log(x$lambda)
           iname <- "negative Log Lambda"
         },
         "dev" = {
           index = x$sumSquare
           iname = "Sum of Squares"
         }
       )
  matplot(index, t(nbeta), lty = 1, xlab = xlab, ylab = ylab,
          type = 'l', ...)
}
