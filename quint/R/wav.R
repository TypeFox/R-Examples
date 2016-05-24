wav <-
function(we,me){
   #  Computes weighted average
   #we =vector of weights
   #me=vector of means
   me<-as.matrix(me)
   index<-!is.na(me)
   wa<-(we[index] %*%me[index])/sum(we[index])
   return(wa) }
