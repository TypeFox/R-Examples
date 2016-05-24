qq.lnorm <-
function(pl,mu,sigma,aveple=TRUE,...){
yvalue <- names(pl)[1]
yvalue<- deparse(substitute(pl))
n<- pl$n[length(pl$n)]
nd<- dim(pl)[1]
ym <- pl[,1]        # data on the vertical axis
  # plotting position
if(aveple){
   pp<-1:nd 
   pp[1]<-pl$ple[1]
   for (j in 2:nd) pp[j]<- (pl$ple[j-1] + pl$ple[j])/2.0
}
else{ pp<- n*pl$ple/(n+1) }

xq<- qnorm(pp)  # normal quantiles on horizontal axis
Rsq <- round(cor(xq,log(ym) )^2,3)

#
plot( xq,ym,type = "n",xlab="Normal Quantile",ylab=yvalue,log="y",... )
points(xq,ym, pch = 1, cex = 0.6,col="red")
#
#  add line

if( missing(mu)|| missing(sigma) ){ fit<- lm( log(ym) ~ xq)
  mu<- fit$coef[1] ; sigma<-fit$coef[2]
   }
GM<- exp(mu  ) ; GSD<- exp(sigma )
xx<- c(xq[1],xq[length(xq)]);  yhat<- exp(mu + xx*sigma)
lines(cbind(xx,yhat),type="b")


par<- c(mu,sigma,Rsq); names(par)<-c("mu","sigma","Rsq")

invisible( list(x=xq,y=ym,pp=pp,par=par) )

}

