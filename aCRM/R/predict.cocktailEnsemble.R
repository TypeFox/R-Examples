predict.cocktailEnsemble <-
function(object,newdata=NULL, ...){
  

predGLM <- predict(object[[1]], newdata, type="response")

predrF <- predict(object[[2]],newdata,type="prob")[,2]

predab <- predict(object[[3]],newdata,type="probs")[,2]

predkF <- predict(object[[4]],newdata,type="probs")


rowMeans(data.frame(predGLM,predrF,predab,predkF))
}
