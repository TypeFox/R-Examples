cpuekapp<-function(x=NULL,nlarge=NULL,absdif=0.001){
   if(is.null(x)) stop ("x is missing")
   if(!is.numeric(x)) stop ("data (x) are not numeric.")
   if(is.null(nlarge)) stop("nlarge must be specified")
   if(!is.numeric(nlarge)) stop ("nlarge must be numeric.")
   ox<-x
   x<-sort(x)
   n<-length(x)
   gu<-function(u){
     a<-0
    for(dd in 1:n) a<-a+dlogis(u,y[dd],hest)/n #checked against equation 3
      a
  }
  Gu<-function(u){
    a<-0
    for(dd in 1:n) a<-a+plogis(u,y[dd],hest)/n # checked against Equation 4
    a
  }
   integrand<-function(u,r){
     u*(Gu(u)^(r-1))*((1-Gu(u))^(n-r))*gu(u)
  }
  E <- function(r) {
   exp((1/beta(r,n-r+1))*integrate(integrand,-Inf,Inf,r)$value)
  }
 
  ssq<-function(h){
     # Calculate maximize
     sgy<-NULL
     sgny<-NULL
      for(i in 1:c(n-nlarge)){
       a<-0
       for(k in 1:n){
         if(k!=i) a<-a+dlogis(y[i],y[k],h)
       }
        sgy[i]<-a/(n-1)
      }      
      a<-0
      for(k in 1:n){
       if(k!=c(n-nlarge)) a<-a+plogis(y[n-nlarge],y[k],h)
      }
      sgny<-a/(n-1)
      sum(log(sgy))+nlarge*log(1-sgny)
   }
   ## Calculations
  m2<-0
  m1<-mean(x)
  while(abs(c(m1-m2))>absdif){
     m1<-mean(x)
     y<-log(x)
     hest<-optimize(ssq,c(1e-10,1000),maximum=TRUE)$maximum 
     for(j in 1:nlarge) x[n-j+1]<-E(n-j+1) 
     m2<-mean(x)
  }  
if(c(m1-m2)<0) outmean<-m1
if(c(m1-m2)>=0) outmean<-m2
     ans<-as.data.frame(matrix(c(c(c(length(x)-nlarge+1):length(x)),ox[c(length(ox)-nlarge+1):c(length(ox))],x[c(length(x)-nlarge+1):c(length(x))]),nrow=nlarge,ncol=3))
     names(ans)<-c("obs","original","expectation")
  output<-list(outmean,ans)
  names(output)<-c("kappmean","expectations")
  return(output)
}# End function
