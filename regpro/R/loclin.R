loclin<-function(arg,x,y,h=1,kernel="gauss",type=0)
{
d<-length(arg)
n<-length(y)

if (d>1){

if (kernel=="bart") 
   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( (2*pi)^(-d/2)*exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ ans<-(rowSums(xx^2) <= 1) 
                      return( ans ) }

argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
w<-ker((x-argu)/h)/h^d
weights<-w/sum(w)

X<-cbind(matrix(1,n,1),x-argu)
W<-diag(weights)
A<-t(X)%*%W%*%X     
invA<-solve(A,diag(rep(1,d+1))) 
B<-t(X)%*%W%*%y
esti<-invA%*%B
est<-esti[type+1]

}
else{  # d==1  #########################################

if (kernel=="gauss") ker<-function(xx){ return( exp(-xx^2/2) ) }
if (kernel=="uniform") ker<-function(xx){ return( (abs(xx) <= 1) ) }

x<-matrix(x,length(x),1)
w<-ker((x-arg)/h)/h^d   
weights<-w/sum(w)

X<-cbind(matrix(1,n,1),x-arg)
W<-diag(c(weights))
A<-t(X)%*%W%*%X     
invA<-solve(A,diag(rep(1,d+1))) 
B<-t(X)%*%W%*%y
esti<-invA%*%B
est<-esti[type+1] 

other<-FALSE
if (other){
w<-ker((arg-x)/h); p<-w/sum(w)
barx<-sum(p*x); bary<-sum(p*y)
q<-p*(1-((x-barx)*(barx-arg))/sum(p*(x-barx)^2))


s1<-sum(w*(x-arg))
s2<-sum(w*(x-arg)^2)
q<-w*(s2-(x-arg)*s1)/sum(w*(s2-(x-arg)*s1))
}

}

return(est)
}

