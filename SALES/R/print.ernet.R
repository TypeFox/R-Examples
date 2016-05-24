print.ernet <- function(x, digits = max(3, getOption("digits") - 
    3), ...) {
    cat("\nCall: ", deparse(x$call), "\n\n")
    print(cbind(Df = x$df, Lambda = signif(x$lambda, digits)))
} 
