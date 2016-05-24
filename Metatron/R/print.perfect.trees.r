print.perfect.trees<-function(x,...){
      if (!is.element("perfect.trees", class(x))) 
        stop("Argument 'x' must be an object of class \"perfect.trees\".")
cat("There are",x$ntrees,"independent trees and", dim(x$coefficients)[1],"parameters in this model. \n")

cat("\n")
cat("Estimation of the prevalences in primary studies and accuracy indices of the test:\n")
print(x$coefficients)
cat("\n")
cat("Model fit statistics:\n")
cat("AIC=",x$aic,"\n")
cat("\n")
print(x$gof)
cat("\n")
}
