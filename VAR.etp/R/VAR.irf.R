VAR.irf <-
function(b,p,sigu,h=10,graphs=FALSE)
{
k <- nrow(b)
mf <- VAR.mainf(b,p,h)
q <- t(chol(sigu))
index <- 1:k
impmat <- matrix(0,nrow=k*k,ncol=h+1)
for( i in 1:(h+1) ){
impmat[,i] <- t( t( as.vector(mf[,index] %*% q)) )
index <- index+k
}
colnames(impmat) <- paste("h",0:h,sep="")
tem1 <- rownames(b)
tem2 <- character()
for (i in 1:length(tem1))
for (j in 1:length(tem1))
tem2 <- c(tem2,paste(tem1[i],tem1[j],sep="->"))
rownames(impmat) <- tem2

impmat=t(impmat)
if (graphs==TRUE){
par(mfrow=c(1,1))
for (i in 1:ncol(impmat)){
plot(impmat[,i],type="l",main=colnames(impmat)[i],ylab="",xlab="h")
abline(h=0,col="red")
par(ask=TRUE) }
}
return(impmat)
}
