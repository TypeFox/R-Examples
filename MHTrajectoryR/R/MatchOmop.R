MatchOmop <- function(signaux, rsStatus){
  who <- which(signaux$AE == "OtherAE")
  if (length(who)>0){
    signaux <- signaux[-who,]
  }
  signaux$statu <- rep(NA, nrow(signaux))
  for (ae in unique(signaux$AE)[1:4]){
    reference <- rsStatus[which(rsStatus$AE==ae),]
    reference <- reference[which(is.na(reference$ATC)==FALSE),]
    who <- which(signaux$AE==ae)
    for (j in who){
      if (any(reference$ATC == signaux$ATC[j])){
        signaux$statu[j] <- reference$statu[which(reference$ATC==signaux$ATC[j])[1]]
      }
      
    }
    mini <- signaux$statu[who]
  }
  return(signaux)
}