seqmon<-function(a,b,t,int){
#pL<-Cumulative probability of y(k) going below a(k)
#pU<-Cumulative probability of y(k) going above b(k)
#y is standard normal under Ho
#y<-z(t)/sqrt(t) where Cov(z(t_1),z(t_2))<-min(t_1,t_2);
#t<-information, usually 1,2,3...k, int<-#intervals for integration
 ones<-function(a,b){array(rep(1,a*b),c(a,b))}
 normcdf<-function(xx){pnorm(xx)}
 d<-(b-a)/int;m<-length(a)
pU=ones(m,1)
pL=ones(m,1) 
sq2pi<-sqrt(2*pi)
H<-1:int[1]
E<-ones(1,int[1])
xo<-a[1]+((1:int[1])-.5*E)*d[1]
pU[1]<-normcdf(-(sqrt(t[1])*b[1])/sqrt(t[1]))
M<-t((d[1]/sq2pi)*exp(-(sqrt(t[1])*xo)^2/(2*t[1])))
pL[1]<-normcdf(sqrt(t[1])*a[1]/sqrt(t[1]))
for (k in 2:m) {
   VU<-normcdf(-(sqrt(t[k])*b[k]*E-sqrt(t[k-1])*xo)/sqrt(t[k]-t[k-1]))
   VL<-normcdf((sqrt(t[k])*a[k]*E-sqrt(t[k-1])*xo)/sqrt(t[k]-t[k-1]))
   pL[k]<-pL[k-1]+VL%*%M
   pU[k]<-pU[k-1]+VU%*%M
   x<-a[k]+((1:int[k])-.5*ones(1,int[k]))*d[k]
   M<-(d[k]*sqrt(t[k])/(sq2pi*sqrt(t[k]-t[k-1])))*
   exp(-(sqrt(t[k])*(t(x)%*%ones(1,int[k-1]))-sqrt(t[k-1])*
   (ones(int[k],1)%*%xo))^2/(2*(t[k]-t[k-1])))%*%M
   xo<-x}
  c(pL,pU)
}


   
