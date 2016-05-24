print.fregre.fd<-function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\n-Call: ", deparse(x$call), "\n", sep = "")
    if (length(coef(x))) {
        cat("\n-Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2,
            quote = FALSE)
      if (x$call[[1]]=="fregre.plm")    cat("\n-Bandwidth (h): ",x$h.opt)
      else if (x$call[[1]]=="fregre.lm")      print(x$beta.est[[2]])
    }
    else {
         if (x$call[[1]]=="fregre.np" || x$call[[1]]=="fregre.np.cv") {
            cat("\n-Bandwidth (h): ",x$h.opt)
                }
         }
    cat("\n-R squared: ",x$r2)
    cat("\n-Residual variance: ",x$sr2,"\n")
    invisible(x)
}

