# Set method to print the output Author: Yiwen Zhang

setMethod("print", signature = "MGLMfit", function(x, ...) {
    print(cbind(estimate = x$estimate, SE = x$SE))
    cat("\n")
    cat("Distribution: ", x$distribution, "\n", sep = "")
    cat("Log-likelihood: ", x$logL, "\n", sep = "")
    cat("BIC: ", x$BIC, "\n", sep = "")
    cat("AIC: ", x$AIC, "\n", sep = "")
    cat("LRT test p value: ", ifelse(x$LRTpvalue > 1e-04, sprintf("%.3f", x$LRTpvalue), 
        "<0.0001"), "\n", sep = "")
    cat("Iterations: ", x$iter, "\n", sep = "")
})



setMethod("print", signature = "MGLMreg", function(x, ...) {
    cat("Call: ")
    print(x$call)
    cat("\n")
    cat("Coefficients:\n")
    print(x$coefficients)
    cat("\n")
    cat("Hypothesis test: \n")
    print(x$test)
    cat("\n")
    cat("Distribution: ", x$distribution, "\n", sep = "")
    cat("Log-likelihood: ", x$logL, "\n", sep = "")
    cat("BIC: ", x$BIC, "\n", sep = "")
    cat("AIC: ", x$AIC, "\n", sep = "")
    cat("Iterations: ", x$iter, "\n", sep = "")
})




setMethod("print", signature = "MGLMsparsereg", function(x, ...) {
    cat("Call: ")
    print(x$call)
    cat("\n")
    cat("Distribution: ", x$distribution, "\n", sep = "")
    cat("Log-likelihood: ", x$logL, "\n", sep = "")
    cat("BIC: ", x$BIC, "\n", sep = "")
    cat("AIC: ", x$AIC, "\n", sep = "")
    cat("Degrees of freedom: ", x$Dof, "\n", sep = "")
    cat("Lambda: ", x$lambda, "\n", sep = "")
    if (!is.null(x$maxlambda)) {
        cat("Max lambda: ", x$maxlambda, "\n", sep = "")
    }
    cat("Iterations: ", x$iter, "\n", sep = "")
})


setMethod("print", signature = "MGLMtune", function(x, ...) {
    cat("Call: ")
    print(x$select$call)
    cat("\n")
    cat("Distribution: ", x$select$distribution, "\n", sep = "")
    cat("Log-likelihood: ", x$select$logL, "\n", sep = "")
    cat("BIC: ", x$select$BIC, "\n", sep = "")
    cat("AIC: ", x$select$AIC, "\n", sep = "")
    cat("Degrees of freedom: ", x$select$Dof, "\n", sep = "")
    cat("Lambda: ", x$select$lambda, "\n", sep = "")
    cat("Number of grid points: ", nrow(x$path), "\n", sep = "")
}) 
