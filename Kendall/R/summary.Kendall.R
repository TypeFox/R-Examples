"summary.Kendall" <-
function(object, ...){
cat(paste(
   "Score = ", object$S,    
    ", Var(Score) =", prettyNum(object$varS)), 
    fill = TRUE)
cat(paste("denominator = ", prettyNum(object$D)), fill=TRUE)
cat(paste(
   "tau = ", format(object$tau, digits = 3),    
    ", 2-sided pvalue =", format.pval( object$sl), sep = ""), 
    fill = TRUE)

}

