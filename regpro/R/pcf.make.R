pcf.make<-function(func,N,support=NULL)
{
d<-length(N)

if (d>1){

recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

if (func=="phi"){

  phi<-function(x){ return( (2*pi)^(-d/2)*exp(-sum(x^2)/2) ) }
 
  if (is.null(support)){
       support<-matrix(0,2*d,1)
       for (i in 1:d){
           support[2*i-1]<--3
           support[2*i]<-3
       }
  }
  lowsuppo<-matrix(0,d,1)
  for (i in 1:d) lowsuppo[i]<-support[2*i-1]
  step<-matrix(0,d,1)
  for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]


  for (i in 1:recnum){
     inde<-digit(i-1,N)+1
     arg<-lowsuppo+step*inde-step/2   # arg<-matrix(arg,d,1)
 
     valli<-phi(arg)    

     value[i]<-valli
     index[i,]<-inde
  }

}


down<-index-1
high<-index

pcf<-list(
value=value,index=index,
down=down,high=high,  
support=support,N=N)

}

return(pcf)
}

