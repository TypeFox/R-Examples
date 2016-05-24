prepareX <-
function(DM_kov){
  
  for(i in 1:ncol(DM_kov)){
    if(is.factor(DM_kov[,i])){
      if(is.ordered(DM_kov[,i])){
        DM_kov[,i] <- as.numeric(DM_kov[,i])
      } else{
        if(nlevels(DM_kov[,i])==2){
          DM_kov[,i] <- as.numeric(DM_kov[,i])-1
        } else{
          labs <- levels(DM_kov[,i])
          DM_kov[,i] <- as.numeric(DM_kov[,i])
          for(j in 2:length(labs)){
            DM_kov[,labs[j]] <- ifelse(DM_kov[,i]==j,1,0)
          }
          DM_kov <- DM_kov[,-i]
        }
      }
    }
  }
  return(DM_kov)
}
