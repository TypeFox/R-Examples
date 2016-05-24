print.cpernet <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall: ", deparse(x$call), "\n\n")
    print(cbind(Df1 = x$df.beta, Df2 = x$df.theta,
                Lambda = signif(x$lambda, digits)))
} 
