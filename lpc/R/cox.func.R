cox.func <- function(x,y,censoring.status,s0=0){
  scor <- coxscor(x,y, censoring.status)$scor
  sd <- sqrt(coxvar(x,y, censoring.status))
  tt <- scor/(sd+s0)
  return(list(tt=tt, numer=scor, sd=sd))
}

