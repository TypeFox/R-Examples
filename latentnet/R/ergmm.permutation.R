ergmm.permutation <- function(n)
{
  if(n==2)
    return(matrix(c(1,2,2,1),2,2))
  temp <- ergmm.permutation(n-1)
  for(i in 0:(n-1))
  {
    if(i==0)
      temp2 <- cbind(temp,n)
    else{
      if (i==(n-1))
        temp2 <- rbind(temp2,cbind(n,temp))
      else
        temp2 <- rbind(temp2,cbind(temp[,1:i],n,temp[,(i+1):(n-1)]))
    }
  }
  colnames(temp2)<-1:n
  return(temp2)
}

which.perm.nearest<-function(K.ref,K){
  perms<-ergmm.permutation(max(c(K.ref,K)))
  order(as.numeric(perms[which.min(apply(perms[,K],1,function(K.perm)sum(K.perm!=K.ref))),]))
}
