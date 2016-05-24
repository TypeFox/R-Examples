VAR.ys <-
function(x,b,p,e,type)
{
n <- nrow(x); k <- nrow(b)

b0 <- b[,1:(k*p),drop=FALSE]
if(type=="const") b1 <- as.matrix(b[,(k*p)+1],nrow=k)
if(type=="const+trend")
{b1 <- as.matrix(b[,(k*p)+1],nrow=k); b2 <- as.matrix(b[,(k*p)+2],nrow=k)}

y <- x[1:p, ,drop=FALSE]

for( i in (p+1):n )
{
index <- 1:k
ytem <- y[nrow(y):1, ,drop=FALSE]
    d1 <- 0
    for(j in 1:p)
    {
    d <- b0[,index] %*% t(ytem[j, , drop=FALSE])
    d1<- d1+d
    index<- index+k
    }

d1 <- d1 + as.matrix(e[i-p,],nrow=k)
if(type=="const") d1 <- d1 + b1
if(type=="const+trend") d1 <- d1 + b1 + b2*i
colnames(d1) <- NULL
y <- rbind(y,t(d1))
}

return(y)
}
