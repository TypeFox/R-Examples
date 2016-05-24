disc.integrate <-
function(f,z0,R=0,gamma=c(0,2*pi),rel.tol= .Machine$double.eps^0.5,subdivisions=100L){
if(R<0 | is.infinite(R)){stop("Invalid radius")}
Result<-NULL

Integral<-function(f0,eval.pt,radius,start,end,error,interval){
U<-function(theta){Re(f0(eval.pt+radius*exp(1i*theta)))*cos(theta)-Im(f0(eval.pt+radius*exp(1i*theta)))*sin(theta)}
V<-function(theta){Re(f0(eval.pt+radius*exp(1i*theta)))*sin(theta)+Im(f0(eval.pt+radius*exp(1i*theta)))*cos(theta)}
Re<-integrate(U,start,end,rel.tol=error,subdivisions=interval)
Im<-integrate(V,start,end,rel.tol=error,subdivisions=interval)
answer<-radius*1i*Re$value-radius*Im$value
abs.error<-sqrt(Re$abs.error^2+Im$abs.error^2)
return(list(value=answer,abs.error=abs.error))}

if(R>0){I<-Integral(f0=f,eval.pt=z0,radius=R,start=gamma[1],end=gamma[2],error=rel.tol,interval=subdivisions)
Result<-list(value=I$value,abs.error=I$abs.error)}

if(R==0){
record<-ratio<-NULL;eps<-1e-12
repeat{
I<-NULL
repeat{
tryCatch(I<-Integral(f0=f,eval.pt=z0,radius=eps,start=gamma[1],end=gamma[2],error=rel.tol,interval=subdivisions), error=function(e) NULL)
if(length(I)>0 & eps==1e-12 & length(record)<1) {eps<-1e-307; record<-NULL; I<-NULL; break}
if(length(I)>0 | eps>1) break
eps<-eps*10;if(eps>1) break}
if(length(I)>0){if(length(record)>0){ratio<-c(ratio,Mod(record[nrow(record),1]/I$value))}else{ratio<-c(ratio,NA)};
record<-rbind(record,c(I$value,I$abs.error,eps))}
eps<-eps*10
if(eps>1) break}
#test stability
if(length(record)>0){
value<-record[1,1];abs.error<-record[1,2]
colnames(record)<-c("Value","Abs.error","Radius")
#if(min(Re(ratio[is.finite(ratio)]))>1) warning("Result may not be reliable (value is sensitive to the radius).
#It might be that the evaluation point of the integral is a branch point")
Result<-list(value=value,abs.error=abs.error,record=record)}}
return(Result)}
