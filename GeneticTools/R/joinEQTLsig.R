joinEQTLsig <- function(eqtlTemp){
  output <- eqtlTemp[[1]]
  if(length(eqtlTemp)>1)
  {
    for(i in 2:length(eqtlTemp)){
      output <- rbind(output,eqtlTemp[[i]])
    }
  }
  output
}