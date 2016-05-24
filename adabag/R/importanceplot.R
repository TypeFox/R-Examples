importanceplot <-
function(object,...){
  if(!((class(object)=="bagging")|(class(object)=="boosting")))
    stop("object class should be bagging or boosting")
  
  barplot(object$imp[order(object$imp,decreasing=TRUE)], main="Variables relative importance", col="lightblue",las=1,xaxs="r",...)
}
