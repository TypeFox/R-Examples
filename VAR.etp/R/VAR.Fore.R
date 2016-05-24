VAR.Fore <-
function(x,b,p,h,type="const")
{
n <- nrow(x); k <- ncol(x); 
b1 <- b[,1:(k*p)]; b2 <- 0; b3 <- 0; tm <- matrix(0,nrow=h)
if(type=="const") b2 <- b[,k*p+1,drop=FALSE]; 
if(type=="const+trend") {b2 <- b[,k*p+1,drop=FALSE]; b3 <- b[,k*p+2,drop=FALSE]; tm <- as.matrix((n+1):(n+h))}
rx <- x[(n-p+1):n, ,drop=FALSE]
pp <- p
for( i in 1:h)
    {index <- 1:k
    f <- b2+b3*tm[i]
    for(j in 1:p)
    {f1 <- b1[,index] %*% t(rx[(pp-j+1),,drop=FALSE])
    f <- f+f1;
    index <- index+k}
    rx <- rbind(rx,t(f))
    pp <- pp+1
    }
fore <- rx[(p+1):(p+h),,drop=F]
colnames(fore) <- colnames(x); rownames(fore) <- paste("h",1:h,sep="")
return(fore)
}
