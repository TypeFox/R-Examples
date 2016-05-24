pcf.loclin<-function(x,y,h,N,type=0,kernel="gauss",support=NULL,
alt=FALSE,alt2=FALSE)
{
# type=0 ,jos regfunc, type=1, jos 1. muuttujan osit.deriv,
d<-length(N)
n<-length(y)

if (d>1){

if (kernel=="bart") 
   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( (2*pi)^(-d/2)*exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ return( (rowSums(xx^2) <= 1) ) } 

if (kernel=="gauss") radi<-2*h else radi<-h

recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

if (is.null(support)){
  support<-matrix(0,2*d,1)
  for (i in 1:d){
     support[2*i-1]<-min(x[,i])
     support[2*i]<-max(x[,i])
  }
}
lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

for (i in 1:recnum){
     inde<-digit(i-1,N)+1
     arg<-lowsuppo+step*inde-step/2

     argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
     w<-ker((x-argu)/h)/h^d
     painot<-w/sum(w)

     X<-cbind(matrix(1,n,1),x-argu)
     W<-diag(painot)
     A<-t(X)%*%W%*%X     
     invA<-solve(A,diag(rep(1,d+1))) 
     B<-t(X)%*%W%*%y
     valli<-invA%*%B

     value[i]<-valli[1+type]
     index[i,]<-inde
}

down<-index-1
high<-index

pcf<-list(
value=value,index=index,
down=down,high=high,  
support=support,N=N)

}
else{  # d==1  #########################################

x<-matrix(x,length(x),1)

if (kernel=="gauss") ker<-function(xx){ return( (2*pi)^(-1/2)*exp(-xx^2/2) ) }
if (kernel=="uniform") ker<-function(xx){ return( (abs(xx) <= 1) ) }
if (kernel=="bart") ker<-function(xx){ (1-xx^2)*(xx^2<=1) }

index<-seq(1:N)
len<-length(index)
down<-matrix(0,len,1)
high<-matrix(0,len,1)
down[,1]<-index-1
high[,1]<-index

value<-matrix(0,N,1)
if (is.null(support)){
   support<-matrix(0,2,1)
   support[1]<-min(x)
   support[2]<-max(x)
}
step<-(support[2]-support[1])/N
lowsuppo<-support[1]

for (i in 1:N){
     inde<-i
     argu<-lowsuppo+step*inde-step/2
     w<-ker((x-argu)/h)/h
     painot<-w/sum(w)
     #if (!is.null(painot)) value[i]<-t(painot)%*%w else value[i]<-mean(w)

     X<-cbind(matrix(1,n,1),x-argu)
     W<-diag(c(painot))   # huom matriisi muutetaan jonoksi
     A<-t(X)%*%W%*%X     
     invA<-solve(A,diag(rep(1,d+1))) 
     B<-t(X)%*%W%*%y
     valli<-invA%*%B

     value[i]<-valli[1+type]
     index[i]<-inde     

     if (alt2==TRUE){
       s2<-c(t(painot)%*%x^2)
       s1<-c(t(painot)%*%x)
        if (type==0){
           q0<-painot*(s2-s1*x)/(s2-s1^2)
           q1<-(painot*x-s1)/(s2-s1^2)
           value[i]<-t(q0)%*%y+t(q1)%*%y*argu
        }
        else{
           q<-(painot*x-s1)/(s2-s1^2)
           value[i]<-t(q)%*%y
        }
         index[i]<-inde   
     }

    if (alt==TRUE){
       s2<-c(t(painot)%*%(x-argu)^2)
       s1<-c(t(painot)%*%(x-argu))
        if (type==0){
           q<-painot*(s2-s1*(x-argu))/(s2-s1^2)
           value[i]<-t(q)%*%y
        }
        else{
           mx<-c(t(painot)%*%x)
           q<-(painot*(x-mx))/(s2-s1^2)
           value[i]<-t(q)%*%y
        }
         index[i]<-inde   
     }

}

pcf<-list(
value=value,
down=down,high=high,
support=support,N=N)

}

return(pcf)
}


