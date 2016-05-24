plie2 <-
function(transmat,par.trans){

  M = dim(transmat)[1]
  Y = matrix(0,length(transmat)+length(par.trans),1)
  for (i in (1:(M-1))){
     Y[(M*(i-1)+1):(M*i)] = transmat[,i]
     }
  Y[(M*(M-1)+1):length(Y)]=c(par.trans)
  return(Y)

}
