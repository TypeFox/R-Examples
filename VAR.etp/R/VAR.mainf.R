VAR.mainf <-
function(b,p,h)
{
k <- nrow(b)
b2 <- b[,1:(k*p),drop=FALSE]
mf <- cbind(diag(k),b2[,1:k])

for( i in 2:h) 
{
index1 <- seq(from=(i-1)*k+1,length.out=k)
index2 <- 1:k
    mf0 <- matrix(0,nrow=k,ncol=k)
    for( j in 1:p) 
    {
    a <- b2[,index2]
        
        if( index1[1] < 0 ) mf1 <- matrix(0,nrow=k,ncol=k) 
        else mf1 <- mf[,index1] %*% a
        
    mf0 <- mf0 + mf1
    index1 <- index1-k 
    index2 <- index2+k
    }
    mf <- cbind(mf,mf0)
}
return(mf)
}
