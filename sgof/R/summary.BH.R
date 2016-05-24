summary.BH <-
function(object, ...){
cat("Call:\n")
 print(object$call)
cat("\n")
cat("Parameters:","\n")
cat("alpha=",object$alpha,"\n")
tabla<-table(object$Adjusted.pvalues<=object$alpha)

if(sum(object$Adjusted.pvalues>object$alpha)==length(object$data)){attributes(tabla)$dimnames[[1]]=c(">alpha")}else{if(sum(object$Adjusted.pvalues<=object$alpha)==length(object$data)){attributes(tabla)$dimnames[[1]]=c("<=alpha")}else{attributes(tabla)$dimnames[[1]]=c(">alpha","<=alpha")}}



cat("\n")
res <- list(Rejections=object$Rejections,FDR=object$FDR,Adjusted.pvalues=tabla)

class(res) <- "summary.BH"
return(res)

}
