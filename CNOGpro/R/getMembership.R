getMembership <-
function(obs,windowlength,HMMtable){
  start <- as.integer(HMMtable$Startpos[1])
  end <- as.integer(HMMtable$Endpos[nrow(HMMtable)])
  window <- seq(start,end,windowlength)
  membership <- numeric(length=length(obs))
  for (row in 1:nrow(HMMtable)){
    condition <- isbetween(window,HMMtable$Startpos[row],HMMtable$Endpos[row])
    membership[condition] <- as.integer(HMMtable$State[row])
  }
  return(membership)
}
