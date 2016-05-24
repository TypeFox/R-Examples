matrix2vector<-function(betam,variant)
{  
  betav<-NULL
  for (k in 1:nrow(betam)){
    if (variant[k]){
      betav<-c(betav,betam[k,])
    } else{
      betav<-c(betav,betam[k,1])
    }
  }
  return(betav)
}

