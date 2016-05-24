ARM2.Fore <-
function(x,y,M,h)
{
a=M$ARc[,1];b=M$IARM[,1];p=nrow(M$ARc)-1

dat = cbind(y,x);
n <- nrow(dat); k <- ncol(dat);

bmat=matrix(0,nrow=k,ncol=p*k+1)
index=2*(1:p)
bmat[1,index]=b[2:(p+1)]; bmat[1,p*k+1]= b[1]
bmat[2,index]=a[1:p]; bmat[2,p*k+1]= a[p+1]
 
b1 <- bmat[,1:(k*p)]; 
b2 <- bmat[,k*p+1,drop=FALSE]; 

rx <- dat[(n-p+1):n, ,drop=FALSE]
pp <- p
for( i in 1:h)
    {index <- 1:k
    f <- b2
    for(j in 1:p)
    {f1 <- b1[,index] %*% t(rx[(pp-j+1),,drop=FALSE])
    f <- f+f1;
    index <- index+k}
    rx <- rbind(rx,t(f))
    pp <- pp+1
    }
fore <- rx[(p+1):(p+h),,drop=F]
colnames(fore) <- c("y","x"); rownames(fore) <- paste("h",1:h,sep="")
return(fore)
}
