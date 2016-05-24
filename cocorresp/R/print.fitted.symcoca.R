"print.fitted.symcoca" <-
function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    cat("\nFitted values for:", x$nam.dat["namY"], "\n")
    print(zapsmall(x$Yhat1), digits = digits, ..., print.gap = 2)
    cat("\nFitted values for:", x$nam.dat["namX"], "\n")
    print(zapsmall(x$Yhat2), digits = digits, ..., print.gap = 2)
    invisible(x)
  }

