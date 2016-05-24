

coef.src <- function(object, ...){
     p3 <- round(t(object$coeff[, "Estimate", drop = FALSE]), 5)
     rownames(p3) <- ""
     print(p3)
     invisible(object)
}