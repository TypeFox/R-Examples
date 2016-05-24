deplie2 <-
function(Y,M=NULL) { 

  if (is.null(M)) {M = floor((-1+sqrt(1+4*length(Y)))/2)}
  trans = matrix(0,M,M-1);
  nc = (length(Y)-(M*(M-1)))/M
  par.trans = matrix(0,M,nc);
  for (i in 1:(M-1)){
     trans[,i] = Re(Y[(M*(i-1)+1):(M*i)]);
     }
  for (i in 1:nc) {
  #	par.trans[,i] = Y[(M*(M-1)+(i-1)*nc+1 ):(M*(M-1)+i*nc )]
  	par.trans[,i] = Y[(M*(M-1)+(i-1)*M+1 ):(M*(M-1)+i*M )]
  }
  res <- NULL
  res$trans <- trans
  res$par <- par.trans
  return(res)

}
