VARB.ys <-
function(x,b,p,e,type)
{
n <- nrow(x); k <- nrow(b)

b0 <- b[,1:(k*p),drop=FALSE]
if(type=="const") b1 <- as.matrix(b[,(k*p)+1],nrow=k)
if(type=="const+trend")
{b1 <- as.matrix(b[,(k*p)+1],nrow=k); b2 <- as.matrix(b[,(k*p)+2],nrow=k)}

y <- x[(n-p+1):n, ,drop=FALSE]
for( i in (n-p):1 )
{
index <- 1:k

    d1 <- 0
    for(j in 1:p)
    {
    d <- b0[,index] %*% t(y[j, , drop=FALSE])
    d1<- d1+d
    index<- index+k
    }
d1 <- d1 + as.matrix(e[i,],nrow=k)
if(type=="const") d1 <- d1 + b1
if(type=="const+trend") d1 <- d1 + b1 + b2*i
colnames(d1) <- NULL
y <- rbind(t(d1),y)
}

return(y)
}
