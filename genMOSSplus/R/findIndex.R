findIndex <-
function (obsVector, dimens) {

  cumDimens <- array(0, length(obsVector))

  cumDimens[1] <- 1
  
  for (i in 2:length(obsVector))
    cumDimens[i] <- cumDimens[i-1] * dimens[i-1];    
 
  index <- sum(obsVector * cumDimens) + 1

  return(index)  

}
