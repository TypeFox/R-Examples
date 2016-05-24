topological.binary.main=function(x,B,alpha=0.05,critical.value){
  n=length(x)
  k=ceiling(n/B)
  
  blocks=array(NA,dim=c(k,B))
  blocks=matrix(unlist(split(x,k)),ncol=B,byrow = TRUE)  
 
  num.different=nrow(unique(blocks))
  if (num.different<1){
    sonucTBT=0
  }else if (num.different>=min(k,2^B)){
    sonucTBT=1
  }else if (num.different<critical.value){
      sonucTBT=0
  }else{
      sonucTBT=1
  }

  result=list(result.TBT=sonucTBT,statistic=num.different,name="TBT")
  return(result)  
}