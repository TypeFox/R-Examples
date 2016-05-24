#############################################################
#
#	print.quantile.localdepth function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: September, 2, 2008
#	Version: 0.1
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

print.quantile.localdepth <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    if (is.list(x)) {
      cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
      cat("Quantile: \n")
      print.default(format(x$quantile, digits=digits),
                  print.gap = 2, quote = FALSE)
      cat("Statistics on the dimension: \n")
      print(summary(x$stats), digits=digits, ...)
    } else {
      print.default(format(x, digits=digits),
                  print.gap = 2, quote = FALSE)
    }
    cat("\n")
    invisible(x)
}
