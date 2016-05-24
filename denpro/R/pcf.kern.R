pcf.kern<-function(dendat,h,N,kernel="gauss",weights=NULL,support=NULL,
lowest=0,radi=0)
{
d<-length(N)

if (d>1){

if (length(h)==1) h<-rep(h,d)

if (kernel=="bart") 
   ker<-function(xx,d){ 
         musd<-2*pi^(d/2)/gamma(d/2)
         c<-d*(d+2)/(2*musd)
         return( c*(1-rowSums(xx^2))*(rowSums(xx^2) <= 1) ) 
   }
if (kernel=="gauss") 
   ker<-function(xx,d){ return( (2*pi)^(-d/2)*exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx,d){ 
         c<-gamma(d/2+1)/pi^(d/2) 
         return( (rowSums(xx^2) <= 1) ) 
   } 
if (kernel=="epane") 
   ker<-function(xx,d){ 
      c<-(3/4)^d 
      xxx<-(1-xx^2)*(1-xx^2>=0)
      return( c*apply(xxx,1,prod) ) 
   } 

if (is.null(radi)) if (kernel=="gauss") radi<-2*h else radi<-h

recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

if (is.null(support)){
  support<-matrix(0,2*d,1)
  for (i in 1:d){
     support[2*i-1]<-min(dendat[,i])-radi
     support[2*i]<-max(dendat[,i])+radi
  }
}
lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

numpositive<-0
for (i in 1:recnum){
     inde<-digit(i-1,N)+1
     arg<-lowsuppo+step*inde-step/2
     argu<-matrix(arg,dim(dendat)[1],d,byrow=TRUE)
#     neigh<-(rowSums((argu-x)^2) <= radi^2)
#     if (sum(neigh)>=2){     # if there are obs in the neigborhood
#
#       xred<-dendat[neigh,]
#       argu<-matrix(arg,dim(xred)[1],d,byrow=TRUE)

       xxx<-sweep(dendat-argu,2,h,"/")
       w<-ker(xxx,d)/prod(h)
       valli<-mean(w)
       if (!is.null(weights)) valli<-t(weights)%*%w
#     }
#     else valli<-mean(y)

      if (valli>lowest){
           numpositive<-numpositive+1
           value[numpositive]<-valli
           index[numpositive,]<-inde
      }
      #value[i]<-valli
      #index[i,]<-inde

}

value<-value[1:numpositive]
index<-index[1:numpositive,]
down<-index-1
high<-index

pcf<-list(
value=value,index=index,
down=down,high=high,  
support=support,N=N)

}
else{  # d==1  #########################################

d<-1
x<-matrix(dendat,length(dendat),1)

if (kernel=="gauss") ker<-function(xx,d){ return( (2*pi)^(-1/2)*exp(-xx^2/2) ) }
if (kernel=="uniform") ker<-function(xx,d){ return( (abs(xx) <= 1) ) }

index<-seq(1:N)
len<-length(index)

value<-matrix(0,N,1)
if (is.null(support)){
   support<-matrix(0,2,1)
   support[1]<-min(x)
   support[2]<-max(x)
}
step<-(support[2]-support[1])/N
lowsuppo<-support[1]

numpositive<-0
for (i in 1:N){
     inde<-i
     argu<-lowsuppo+step*inde-step/2
     w<-ker((x-argu)/h,1)/h
     if (!is.null(weights)) valli<-t(weights)%*%w else valli<-mean(w)
     if (valli>lowest){
           numpositive<-numpositive+1
           value[numpositive]<-valli
           index[numpositive]<-inde
     }
}

value<-value[1:numpositive]
index<-index[1:numpositive]

down<-matrix(0,numpositive,1)
high<-matrix(0,numpositive,1)
down[,1]<-index-1
high[,1]<-index

pcf<-list(
value=value,
down=down,high=high,
support=support,N=N)

}

return(pcf)
}
