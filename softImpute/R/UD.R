UD=function(U,D,n=nrow(U)){
  U*outer(rep(1,n),D,"*")
}
