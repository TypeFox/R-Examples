"print.summary.fitCoinertia" <-
function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    cat("\nCoinertia analysis of the triplets (Q_1,K_1,R_0) and (Q_2,K_2,R_0)\n\n")
    writeLines(strwrap(pasteCall(x$call)))
    cat("\nEigenvalues:\n")
    print(round(x$lambda, digits = digits), ..., print.gap = 3)
    cat("\nU-matrix 1:\n")
    print(x$U1, digits = digits, ..., print.gap = 2)
    cat("\nU-matrix 2:\n")
    print(x$U2, digits = digits, ..., print.gap = 2)
    cat("\nX-matrix 1:\n")
    print(x$X1, digits = digits, ..., print.gap = 2)
    cat("\nX-matrix 2:\n")
    print(x$X2, digits = digits, ..., print.gap = 2)
    invisible(x)
  }

