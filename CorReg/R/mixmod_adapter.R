# 'adapting mixmod object to be more efficient in C++
# ' @param list from calcul_mixmod 
mixmod_adapter<-function(mixmod){
  res=mixmod
  if(!is.null(res$details)){
    mat=res$details[[1]]
    if(length(res$details)>1){
      for (i in 2:length(res$details)){
        mat=rbind(mat,res$details[[i]])
      }
    }
    res$details=mat
  }
  return(res)
}