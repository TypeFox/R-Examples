RazonVero <-
function(lambda,x,n2,a,Z=0){

##########################################
PVerosimilitud=function(lambda,x,n2,a){
ad=sum(x); alpha=lambda/a[2]; bet=-a[1]/a[2]

c0=alpha*(1.0-alpha)*x[1];
c1=bet*ad-n2[1]*alpha*(1.0-alpha)-alpha*bet*(n2[2]+2*x[1])
c2=bet*(alpha*(n2[2]+2*n2[1])-(n2[1]+x[2])-bet*(n2[2]+x[1]) );
c3=bet^2*sum(n2)

A=-4.5*c3*(c1*c2-3.0*c0*c3)+c2*c2*c2; B=c2*c2-3.0*c1*c3;
u=sqrt(B);ii=A<0; u[ii]=-u[ii]
w1= A/(u^3); w1[w1>1]=1; w1[w1<(-1)]=-1
fhi=(pi+ acos(w1))/3.0
	
pvm1=(-1*c2+2*cos(fhi)*u)/(3*c3); pvm1[pvm1>1]=1;pvm1[pvm1<0]=0
pvm2=alpha+bet*pvm1;pvm2[pvm2>1]=1;pvm2[pvm2<0]=0

cbind(pvm1,pvm2)
}



pvm=PVerosimilitud(lambda,x=x,n=n2,a=a)
pm=x/n2
gk1=log(pvm[,1]/pm[1]); gk2=log((1-pvm[,1])/(1-pm[1]))
gk3=log(pvm[,2]/pm[2]); gk4=log((1-pvm[,2])/(1-pm[2]))
aux=-2*(x[1]*gk1+(n2[1]-x[1])*gk2+x[2]*gk3+(n2[2]-x[2])*gk4)
as.numeric(aux-Z^2)}
