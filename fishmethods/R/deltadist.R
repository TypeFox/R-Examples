 deltadist <- function(x=NULL) {
    if(is.null(x)) stop ("catch data are missing.")
    if(!is.numeric(x)) stop ("catch data are not numeric.")
       G<-function(t){
            j<-2
            s<-(((m-1)**3)/((m+1)*(m**2)*2))*t**2
            a<-1+(((m-1)/m)*t)+s
            while (s>0.000001){
            j<-j+1
            b<-(((m-1)**2)/(m*(m+2*j-3)*j))*t
            s<-s*b
            a<-a+s
            }
            a
           }
   x<-x[!is.na(x)]
   n<-length(x)
   m<-length(x[x>0])
  
if(m>1){
     my<-mean(log(x[x>0]))
     vary<-var(log(x[x>0]))
     c<-(m/n)*exp(my)*G(0.5*vary)
     varc<-(m/n)*exp(2*my)*((m/n)*(G(0.5*vary)^2)-(((m-1)/(n-1))*G(((m-2)/(m-1))*vary)))
     }
if(m==1){
    c<-x[x>0]/n
    varc<-(x[x>0]^2)/n
  }
if(m==0){
    c<-0
    varc<-0
  }
outpt<-unlist(list(deltamean=c,var.mean=varc))
 return(outpt)
}
