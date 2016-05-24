print.modprobs <-
function(x,digits = max(3, getOption("digits") - 3),...){

cat("Posterior model probabilities:\n")
print(x$table,right=FALSE,digits=digits)
cat("\n")
cat("Total number of models visited = ",round(x$totmodsvisit,digits=digits),"\n")}
