"print.fitCoinertia" <-
function(x, axes = c(1:min(6, x$n.axes)),
         digits = max(3, getOption("digits") - 3), ...)
  {
    cat("\nCoinertia analysis of the triplets (Q_1,K_1,R_0) and (Q_2,K_2,R_0)\n\n")
    writeLines(strwrap(pasteCall(x$call)))
    cat("\nEigenvalues:\n")
    print(round(x$lambda, digits = digits), ..., print.gap = 3)
    invisible(x)
  }
