Stine <-
function(bb,n,p)
{
a <- c(1,bb[1:p])
{
if (p == 1)
{
b <- rbind(c(0,0),c(1,3))
}
else
{
value <- 0.5*p - as.integer(0.5*p)

b1 <- matrix(0,nrow=p+1,ncol=p+1)
diag(b1) <- 0:p
b21 <- matrix(0,nrow=p+1,ncol=1)
b22 <- matrix(0,nrow=p+1,ncol=1)
b2 <-  matrix(0,nrow=p+1,ncol=1)

  if (value == 0)
  {
     for (i in 0:(0.5*p-1))
     {
     e <- matrix(0,nrow=p+1,ncol=1)
     index <- seq(i+3,length.out=(0.5*p-i),by = 2)
     e[index] <- matrix(1,nrow=(0.5*p-i),ncol=1)
     b21 <- cbind(b21,-e)
     b22 <- cbind(e,b22)
     }
     b2 <- cbind(b21[,2:(0.5*p+1)],b2,b22[,1:(0.5*p)])
  }  
   else
  {
     for (i in 0:(0.5*(p-1)) ) 
     {
     d <- matrix(0,nrow=p+1,ncol=1) 
     index <- seq(i+2,length.out=( 0.5*(p-1)+1-i ),by = 2)
     d[index] <- matrix(1,nrow=( 0.5*(p-1)+1-i ),ncol=1)
     b21 <- cbind(b21,-d)
     b22 <- cbind(d,b22)
     }
     b2 <- cbind(b21[,3:(0.5*(p-1)+2)],b2,b22[,1:(0.5*(p-1)+1)])
  } 


  b3 <- matrix(0,nrow=p+1,ncol=p+1)
  for (i in 1:(p+1))
  {
    for (j in 1:(p+1))
    {
    if (j < i & i <= p-j+2)
    b3[i,j] <- -1
    if (p-j+2< i & i <= j)
    b3[i,j] <- 1
    }
  }
  
b <- b1+b2+b3
b[,1] <- -b[,1]

}
}

bmat <- -b/n;
bmat1 <- bmat[2:(p+1),1]
bmat2 <- bmat[2:(p+1),2:(p+1)]

eye <- matrix(0,nrow=p,ncol=p)
diag(eye) <- rep(1,p)
ahat <- solve(eye+bmat2) %*% (a[2:(p+1)]-bmat1)
return(ahat)
}
