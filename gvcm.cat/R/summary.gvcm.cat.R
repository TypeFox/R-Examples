summary.gvcm.cat <-
function (
object, dispersion = NULL, ...
)
{
# check input
if (("gvcm.cat" %in% is(object))==FALSE )
     stop ("object must be a 'gvcm.cat' object. \n")

# disp wie in glm
    est.disp <- FALSE
    df.r <- object$df.residual
    if (is.null(dispersion))
        dispersion <- if (object$family$family %in% c("poisson",
            "binomial"))
            1
        else if (df.r > 0) {
            est.disp <- TRUE
            if (any(object$weights == 0))
                warning("observations with zero weight not used for calculating dispersion")
            sum((object$weights * object$residuals^2)[object$weights >
                0])/df.r
        }
        else {
            est.disp <- TRUE
            NaN
        }
       
# index.reduced
index.reduced <- c()
if(dim(object$x.reduction)[2]>0) {
for (i in 1:dim(object$x.reduction)[2]){
   index.reduced <- c(index.reduced,min(which(object$x.reduction[,i]==1))) 
   }
} 

# coefficient table
coef.table <- data.frame(object$coefficients)
coef.table[,1]<- object$coefficients.oml
coef.table[,2]<- object$coefficients
coef.table[index.reduced,3]<- object$coefficients.reduced
coef.table[index.reduced,4]<- object$coefficients.refitted
colnames(coef.table) <- c("coefficients.oml", "coefficients", "coefficients.reduced", 
   "coefficients.refitted" )
   
# df.f?!?   

# defintions
# dev.res <- summary(round(object$residuals,3))
dev.res <- matrix(summary(round(object$residuals,3))[c(1,2,3,5,6)], nrow=1)
colnames(dev.res) <- c("Min", "1Q", "Median", "3Q", "Max")
rownames(dev.res) <- ""

# summary
cat("\nCall:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
    "\n\n", sep = "")

cat("Deviance Residuals: \n")
#cat("    Min        1Q    Median        3Q       Max \n")
#cat(dev.res[1],dev.res[2],dev.res[3],dev.res[5],
#    dev.res[6],"\n \n")
print(dev.res)
cat("\n")

cat("Coefficients: \n")
print(coef.table)
cat("\n")

cat("(Dispersion parameter for ", object$family$family," family taken to be ", 
    dispersion, ") \n", sep="")

cat("    Null deviance: ", object$null.deviance," on ", object$df.null,
    " degrees of freedom \n", sep="")
cat("Residual deviance: ", object$deviance," on ", round(object$df.residual, 2),
    " degrees of freedom \n \n", sep="")

cat("Removed parameters: ", object$number.removed.parameters, " out of ", 
    object$number.selectable.parameters, "\n", sep="")
    
if(object$method %in% c("AIC", "BIC")){
if(object$method %in% c("AIC")){
cat("AIC of chosen model: ", object$tuning, "\n", sep="")
}
if(object$method %in% c("BIC")){
cat("BIC of chosen model: ", object$tuning, "\n", sep="")
}
} else {
#if(object$method %in% c("nlm","lqa")){
cat("Penalization parameter lambda = ", object$tuning[[1]], "\n", sep="")
cat("Tuning: ", "adapted.weights = ", 
    object$control$adapted.weights, ", assured.intercept = ", 
    object$control$assured.intercept, "\n", sep="")
}
cat("Number of iterations: ", object$iter, "\n", sep="")
if (object$converged==TRUE) {cat("The model converged. \n", sep="")} else 
   {cat("The model did not converge. \n", sep="")}

}

  
    
