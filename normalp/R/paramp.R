paramp<-function(x,p)UseMethod("paramp")

paramp.default<-function(x,p=NULL){
if(!is.numeric(x) || (!is.null(p) && !is.numeric(p))) stop (" Non-numeric argument to mathematical function")
ff<-function(Mp) sum(abs(x-Mp)^p)
  Mp<-mean(x)
if (is.null(p)){  
df<-2  
pp<-2
  iter<-0
  i<-0
  p<-estimatep(x,Mp,pp)
  while (abs(p-pp)>0.0001) {
  pp<-p
  op<-optim(Mp,ff,method="BFGS")
  Mp<-op$par
  p<-estimatep(x,Mp,pp)
  i<-i+1
  if (i==100) {iter<-1; break}
  }
  }
else {
if (p<1) stop("Value of p must be greater or equal then 1")
 df<-1  
  op<-optim(Mp,ff,method="BFGS")
  Mp<-op$par
  iter<-0
}
  if (p==1) Mp<-median(x)
  Sp<-((sum(abs(x-Mp)^p))/(length(x)-df))^(1/p)
  if (p>=11.5) {Mp<-mean(c(max(x),min(x)));   Sp<-(max(x)-min(x))/2 }
sd<-sqrt((length(x)-1)*var(x)/length(x))
mn<-mean(x)
ris<-list()
ris$mean<-mn
ris$mp<-Mp
ris$sd<-sd
ris$sp<-Sp
ris$p<-p
ris$iter<-iter
class(ris)<-"paramp"
ris
} 
 
print.paramp<-function(x,...){
dat<-c(Mean=x$mean,Mp=x$mp,Sd=x$sd,Sp=x$sp,p=x$p)
print(dat)
cat("\nno.conv =", as.logical(x$iter),"\n\n")
invisible(x)
}

