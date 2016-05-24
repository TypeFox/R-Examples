## Notation
##    subject specific intervals
##                number: N
##         running index: i
##    support (Peto) intervals
##                number: M
##         running index: m
PetoInt<-function(L,R,status){ 
  #Status: 0 right censored, 1 exact time, 2 interval cencored.
  #R[status==0] <- max(R)+1 #to ensure a right endpoint.
  #it is outcomented because this is done in compGMLE...R instead.
  
  names(L)[status!=1] <- 'L'
  names(R)[status!=1] <- 'R'
  names(L)[status==1] <- 'EL'
  names(R)[status==1] <- 'ER'
  peto.intervals <- c(L,R)
  level.int <- factor(names(peto.intervals),levels=c('R','EL','ER','L'))
  right.order <- order(peto.intervals,level.int)
  peto.intervals <- peto.intervals[right.order]
  tmp1 <- as.numeric(factor(names(peto.intervals), levels=c('R','EL','ER','L')))
  int <- grep('^-3$', diff(tmp1)) #finds the intervals
  tmp2 <- as.numeric(factor(names(peto.intervals), levels=c('EL','R','L','ER')))
  exa <- grep('^3$', diff(tmp2)) #finds the exact observations
  obs.no <- c(int,exa)
  tmp <- peto.intervals[sort(c(obs.no,obs.no+1))]
  out  <-  matrix(tmp,nrow=2)
  out
}
