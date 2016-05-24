CohensD.rawData <-
function(E.data,C.data){  
  d<-(mean(E.data)-mean(C.data))/sqrt(var(C.data)) 
  names(d)<-"d"
  return(d)
}
