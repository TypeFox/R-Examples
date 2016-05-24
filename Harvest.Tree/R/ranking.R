ranking = function(data,numeric.info){
  if(length(numeric.info)>=2){
  order= apply(data[,numeric.info+1],2,rank)
  }
  else if(length(numeric.info==1)){
    order=rank(data[,numeric.info+1])
  }
  data[,numeric.info+1]=order
  
  return(list(data=data,order=order))
    
}