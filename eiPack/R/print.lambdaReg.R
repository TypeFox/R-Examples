print.lambdaReg <- function(x, digits = max(2, getOption("digits") -
                                 4), ...) {
  denom <- as.character(x$call$columns)[2]
  for(i in 3:length(as.character(x$call$columns))){
    denom <- paste(denom,as.character(x$call$columns)[i], sep=", ")
  }
  cat("\nDenominator is sum(", denom,").", "\n\n", sep="")
  print.default(format(x$coef, digits = digits), 
                print.gap = 2, quote = FALSE, ...)
}
