bt.log<-function(meanlog=NULL,sdlog=NULL,n=NULL,alpha=0.05){
    if(is.null(meanlog)|is.null(sdlog)|is.null(n)) stop ("meanlog,sdlog or n are missing")
    if(!is.numeric(meanlog)|!is.numeric(sdlog)|!is.numeric(n)) stop ("meanlog, sdlog, or n are not numeric.")
    sdlog<-sdlog^2   
    G<-function(t){
            j<-2
            s<-(((n-1)**3)/((n+1)*(n**2)*2))*t**2
            a<-1+(((n-1)/n)*t)+s
            while (s>0.000001){
            j<-j+1
            b<-(((n-1)**2)/(n*(n+2*j-3)*j))*t
            s<-s*b
            a<-a+s
            }
            a
         }
     var1<-exp(2*meanlog)*(G(2*sdlog)-G(((n-2)/(n-1))*sdlog))
     outpt<-unlist(list(btmean=exp(meanlog)*G(0.5*sdlog),
     approx.bt.mean=exp(meanlog+sdlog/2),var=var1,sd=sqrt(var1),
     var.mean=var1/n,sd.mean=sqrt(var1/n),
     median=exp(meanlog),LCI=exp(meanlog+sdlog/2-sqrt((sdlog/n)+(sdlog^2)/(2*(n-1)))*qt(1-alpha/2,n-1)),
     UCI=exp(meanlog+sdlog/2+sqrt((sdlog/n)+(sdlog^2)/(2*(n-1)))*qt(1-alpha/2,n-1))))
    return(outpt)
 }

