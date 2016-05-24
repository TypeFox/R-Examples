ARM3.Fore <-
function(x,y,M,h)
{

dat = cbind(y,x);k <- ncol(dat); n <- nrow(dat);
a=M$ARc[,1:(k-1)];b=M$IARM[,1];p=nrow(M$ARc)-1
bmat=matrix(0,nrow=k,ncol=p*k+1)

index1=seq(2,by=k,length.out=p); index2=2:(p+1)
for(i in 1:(k-1)){
bmat[1,index1]=b[index2]; 
index1=index1+1; index2=index2+p}

index1=seq(2,by=k,length.out=p); index2=1
for(i in 2:k){
bmat[i,index1]=a[1:p,index2]; 
index1=index1+1; index2=index2+1}


bmat[1,p*k+1]= b[1]
bmat[2:k,p*k+1]= a[p+1,]

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
forecast <- rx[(p+1):(p+h), ,drop=F]
colnames(forecast) <- c("y",paste("x",1:(k-1),sep="")); rownames(forecast) <- paste("h",1:h,sep="")
return(forecast)
}
