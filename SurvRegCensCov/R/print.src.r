print.src <- function(x, ...){
     
     p1 <- "\nWeibull regression for a right-censored response with an interval-censored covariate\n\n"
     p2 <- "\nCoefficients:\n"
     p3 <- round(t(x$coeff[, "Estimate", drop = FALSE]), 5)
     rownames(p3) <- ""
     p4 <- paste("\nAIC: ", round(x$AIC, 1), "\n", sep = "")
     
     cat(p1)
     cat("Call:\n")
     print(x$Call)
     cat(p2)
     print(p3)
     cat(p4)
     invisible(x)
}