print.gvcm.cat <-
function(x, ...){
cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
    "\n\n", sep = "")
cat("    Null deviance: ", x$null.deviance," on ", x$df.null,
    " degrees of freedom \n", sep="")
cat("Residual deviance: ", x$deviance," on ", round(x$df.residual, 2),
    " degrees of freedom \n \n", sep="")
cat("Removed parameters: ", x$number.removed.parameters, " out of ", 
    x$number.selectable.parameters, "\n", sep="")
if(x$method %in% c("nlm","lqa")){
cat("Penalization parameter lambda = ", x$tuning[[1]], "\n", sep="")
cat("Tuning: ",  "adapted.weights = ", 
    x$control$adapted.weights, ", assured.intercept = ", 
    x$control$assured.intercept, "\n", sep="")
}
if(x$method %in% c("AIC")){
cat("AIC of chosen model: ", x$tuning, "\n", sep="")
}
if(x$method %in% c("BIC")){
cat("BIC of chosen model: ", x$tuning, "\n", sep="")
}

}

