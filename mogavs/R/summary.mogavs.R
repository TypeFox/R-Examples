summary.mogavs <-
function(object, ...){
  ans<-object[c("maxGenerations")]
  
  ans$boundary<-cbind(object$numOfVariables,object$MSE)
  ans$modelsTried<-nrow(object$archiveSet)
  class(ans)<-"summary.mogavs"
  return(ans)
}
