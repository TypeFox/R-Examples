summary.src <- function(object, ...){
     p1 <- "\nWeibull regression for a right-censored response with an interval-censored covariate\n\n"
     p2 <- "\nCoefficients:\n"
     p3 <- paste("\nAIC: ", round(object$AIC, 1), "\n", sep = "")
     cat(p1)
     cat("Call:\n")
     print(object$Call)
     cat(p2)
     options(warn = -1)
     
     c1 <- apply(object$coeff, 1:2, as.numeric)
     c2 <- c1
     c2[, c(1:4, 6:8)] <- format(c2[, c(1:4, 6:8)], digits = 4)
     c2[, 5] <- format.pval(c1[, 5], digits = 3)
     c2 <- data.frame(c2)
     colnames(c2) <- colnames(object$coeff)
     
     print(c2[, 1:5])
     options(warn = 0)
     cat(p3)
     
     class(object) <- "summary.src"
     invisible(object)
}