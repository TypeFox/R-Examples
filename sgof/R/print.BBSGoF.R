print.BBSGoF <-
function(x, ...){
cat("Call:\n") 
print(x$call) 
cat("\n")
cat("Parameters:","\n")
cat("alpha=",x$alpha,"\n")
cat("gamma=",x$gamma,"\n")
cat("kmin=",x$kmin,"\n")
cat("kmax=",x$kmax,"\n")
cat("\n")


if(length(x$deleted.blocks)>1){
cat("Warning:","\n")
cat("Blocks", sort(x$deleted.blocks), "have been removed because they provided negative variances.","\n")
}else if(length(x$deleted.blocks)==1&is.na(x$deleted.blocks)==F&is.infinite(x$deleted.blocks)==F){
cat("Warning:","\n")
cat("Block", sort(x$deleted.blocks), "has been removed because it provided negative variance.","\n")
}

cat("\n")

cat("Rejections:\n") 
print(x$Rejections)
}
