homo.panel <- function(x, panel, window=0.49){
  newxxx <- numeric()
  for(i in 1:length(x$wei)){
    differences <- abs(x$wei[i] - panel)
    v <- which(differences < window)[1]
    if(length(v) > 0){
      y <- panel[v]
    }else{y <- 0}
    newxxx[i] <- y
  }
  x$wei <- newxxx
  return(x)
}