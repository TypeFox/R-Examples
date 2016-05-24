computez.lambdaest <-
function(lambda,nu,max){
  
  forans <- matrix(0,ncol=max+1,nrow=length(lambda))
  for (j in 1:max){
    temp <- matrix(0,ncol=j,nrow=length(lambda))
    for (i in 1:j){temp[,i] <- lambda/(i^nu)}
    for (k in 1:length(lambda)){forans[k,j+1] <- prod(temp[k,])} 
  }
  forans[,1] <- rep(1,length(lambda))
  
  ans <- rowSums(forans)
  
  return(ans)
  
}
