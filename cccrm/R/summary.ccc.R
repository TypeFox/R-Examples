summary.ccc<-
function(object,...){

print(object$model)
cat("\n")
cat("CCC estimated by variance compoments \n")
print(object$ccc[1:4])
}


