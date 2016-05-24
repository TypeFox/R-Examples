vector2matrix<-function(betav,variant,J){  
  betam<-NULL
  index<-1
  for (k in 1:length(variant)){
    if (variant[k]){
      betatemp<-betav[index:(index+J-1)]
      index<-index+J
    } else{
      betatemp<-rep(betav[index],J)
      index<-index+1
    }
    betam<-rbind(betam,betatemp)
  }
  return(betam)
}
