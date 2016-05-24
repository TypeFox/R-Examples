HedgesG.rawData <-
function(E.data,C.data){  
  g<-(mean(E.data)-mean(C.data))/
  sqrt(((length(E.data)-1)*var(E.data)+
  (length(C.data)-1)*var(C.data))/
  (length(E.data)+length(C.data)-2))
  names(g)<-"g"
  return(g)
}
