`BuildMatrix1` <-
function(nodevalues,n){
  # nodevalues = values of the functions f at the nodes
  #              calculated by 'fevaluation'
  # n = number of remaining observed nodes; CAUTION: This 'n' is rather 'nr'
  F <- nodevalues
  s <- length(F[,1])
  N <- matrix(0,nrow=s,ncol=n)
  A0 <- cbind(F,N)
  A0
}

