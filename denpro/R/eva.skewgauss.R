eva.skewgauss<-function(x,mu,sig,alpha)
{

norvak<-prod(sig)^(-1)
point<-(x-mu)/sig
en<-evanor(point)     #dnorm(poi)

point2<-alpha%*%((x-mu)/sig)
pn<-pnorm(point2)

ans<-2*en*pn

return(ans)
}
