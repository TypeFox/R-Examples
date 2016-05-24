single.index.itera<-function(x,y,h=1,M=2,kernel="gauss",vect=FALSE)
{
d<-dim(x)[2]
n<-length(y)
theta0<-matrix(1,d,1)
theta0<-theta0/sqrt(sum(theta0^2))

w<-matrix(0,n,1)
haty<-matrix(0,n,1)
for (m in 1:M){

   xcur<-x%*%theta0
   for (i in 1:n) w[i]<-kernesti.der(xcur[i],xcur,y,h=h,vect=vect)
   weights<-w^2
   for (i in 1:n) haty[i]<-kernesti.regr(xcur[i],xcur,y,h=h,vect=vect)
   z<-xcur+(y-haty)/w

   W<-diag(c(weights),nrow=n,ncol=n)
   A<-t(x)%*%W%*%x     
   invA<-solve(A,diag(rep(1,d))) 
   B<-t(x)%*%W%*%z
   theta0<-invA%*%B
   theta0<-theta0/sqrt(sum(theta0^2))
}

return(theta0)
}

