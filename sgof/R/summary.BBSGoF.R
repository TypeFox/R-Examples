summary.BBSGoF <-
function(object, ...){
cat("Call:\n")
print(object$call)
cat("\n")
cat("Parameters:","\n")
cat("alpha=",object$alpha,"\n")
cat("gamma=",object$gamma,"\n")
cat("kmin=",object$kmin,"\n")
cat("kmax=",object$kmax,"\n")
cat("\n")

if(object$adjusted.pvalues==T){
tabla<-table(object$Adjusted.pvalues<=object$gamma)
if(sum(object$Adjusted.pvalues>object$gamma)==length(object$data)){attributes(tabla)$dimnames[[1]]=c(">gamma")}else{if(sum(object$Adjusted.pvalues<=object$gamma)==length(object$data)){attributes(tabla)$dimnames[[1]]=c("<=gamma")}else{attributes(tabla)$dimnames[[1]]=c(">gamma","<=gamma")}}

res <- list(Rejections=object$Rejections,FDR=object$FDR,Adjusted.pvalues=tabla,Tarone.pvalue.auto=object$Tarone.pvalue.auto, beta.parameters=object$beta.parameters,betabinomial.parameters=object$betabinomial.parameters,sd.betabinomial.parameters=object$sd.betabinomial.parameters,automatic.blocks=object$automatic.blocks)
}else{
res <- list(Rejections=object$Rejections,FDR=object$FDR,Tarone.pvalue.auto=object$Tarone.pvalue.auto, beta.parameters=object$beta.parameters,betabinomial.parameters=object$betabinomial.parameters,sd.betabinomial.parameters=object$sd.betabinomial.parameters,automatic.blocks=object$automatic.blocks)
}
class(res) <- "summary.BBSGoF"
return(res)
}
