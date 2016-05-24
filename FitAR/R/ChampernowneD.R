`ChampernowneD` <-
function(z, p, MeanZero=FALSE){
n<-length(z)
if(MeanZero) y<-z
else y<-z-mean(z)
x0<-x<-y
for (i in 1:p)
    x<-c(x,c(rep(0,i),y[1:(n-i)]))
x<-matrix(x, nrow=p+1, ncol=length(x0), byrow=TRUE)
C<-c(x%*%x0)
A<-toeplitz(C)
E<-matrix(0, nrow=p+1, ncol=p+1)
for (j in 1:p)
    for (i in 1:j){
         E[i+1,j+1] <- E[i,j]+y[i]*y[j]+y[n+1-i]*y[n+1-j]
         }
for (j in 1:(p+1))
    for (i in 1:j)
        E[j,i]=E[i,j]
A-E
}

