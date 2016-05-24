plie2.VM <-
function(transmat,par.trans){

  M = dim(transmat)[1]
  Y = matrix(0,M*(M-1)+2*length(par.trans),1)
  for (i in (1:(M-1))){
     Y[(M*(i-1)+1):(M*i)] = transmat[,i]
     }
  Y[(M*(M-1)+1):(M*(M-1)+length(par.trans))]=c(Re(par.trans))
  Y[(M*(M-1)+length(par.trans)+1):length(Y)]=c(Im(par.trans))
  return(Y)

}
