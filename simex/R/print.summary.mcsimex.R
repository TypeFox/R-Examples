print.summary.mcsimex <-
function (x, digits = max(3, getOption("digits") - 3), ...)
{
cat("Call:\n")
print(x$call)
cat("\nNaive model: \n")
print(x$naive.model)
cat("\nSimex variable :", x$SIMEXvariable, "\n")
cat("Misclassification matrix: \n")
if(is.character(x$mc.matrix))
print(x$mc.matrix) else lapply(x$mc.matrix, print)
cat("\nNumber of iterations: ", x$B, "\n")
cat("\nResiduals: \n")
print(summary(x$residuals), digits)
cat("\nCoefficients: \n")
if (any(names(x$coefficients) == "asymptotic")) {
cat("\nAsymptotic variance: \n")
printCoefmat(x$coefficients$asymptotic, digits = digits)
}
if (any(names(x$coefficients) == "jackknife")) {
cat("\nJackknife variance: \n")
printCoefmat(x$coefficients$jackknife, digits = digits)
}
return(invisible(x))
}

