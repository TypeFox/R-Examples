prune.rt <-
function(tree,data,c.par=NULL,...){
  # tree = rt object
  model<-tree
  # data = dataset on which the full regression trunk model has been estimated
  # c.par = parameter "c" for the "c SE" rule. The value depends on the sample size
  
  if(is.null(summary(model)$REcv)){
    stop("Because estimates of the cross-validated error are lacking, pruning is not possible. Grow the regression trunk again with stima, using the cross-validation procedure.","\n")}
  if(is.null(c.par)){
    if(as.numeric(model$trunk$n[1])<= 100) c.par<-1
    if(as.numeric(model$trunk$n[1])> 100 & as.numeric(model$trunk$n[1])<= 300) c.par<-.8
    if(as.numeric(model$trunk$n[1])> 300) c.par<-.5
  
  }
  cat("The c parameter used in the pruning equals ", c.par,"\n")

  minrow <- which(model$goffull$REcv == min(model$goffull$REcv))[1]
  bestrow <- min(which(model$goffull$REcv <= (model$goffull$REcv[minrow] +c.par * model$goffull$SEcv[minrow])))
   
  if(bestrow==1) {
    cat("The pruned trunk has zero splits --> no interaction terms are present", "\n")
  }
  if(bestrow==2) {
    cat("The pruned trunk has one split --> no interaction terms are present", "\n")
    pruned.trunk<-stima(data,1,first=which(colnames(data)==model$trunk[2,1]),vfold=0)
    pruned.trunk$goffull<-model$goffull[1:bestrow,]
    return(pruned.trunk)
  }
  if(bestrow>2){ 
    cat("The pruned trunk has ", bestrow-1, " splits",  "\n")
           cat("The Cross-Validated Residual Error is", model$goffull[bestrow,5], "\n")
           cat("The Standard Error of the Cross-Validated Residual Error is", model$goffull[bestrow,6], "\n")
           cat("The First Splitting Predictor is: ", model$trunk[2,1], "\n")
           cat("It corresponds to column", which(colnames(data)==model$trunk[2,1]), "\n")
   
      pruned.trunk<-stima(data,bestrow-1,first=which(colnames(data)==model$trunk[2,1]),vfold=0)
      pruned.trunk$goffull<-model$goffull[1:bestrow,]
      return(pruned.trunk)
  }
}
