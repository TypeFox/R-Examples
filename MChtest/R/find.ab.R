"find.ab" <-
function(n=100000,ALPHA=.05,BETA=.2,higha=100){
	
	
higha<-100	
### the function H comes from a one-sided normal test 
### with 80 percent power. This is the cdf
H<-function(x,alpha=.05,beta=.2){
	1-pnorm(
	   qnorm(1-x)-qnorm(1-alpha) - qnorm(1-beta) )
}
### To find the mean of H, we use numeric integration
### Let x={ 0=x_0<x_1<...<x_k=1}
### then 
### \sum_{i=1}^{k-1}  x_i     (H(x_{i+1}) - H(x_i))
###   <=   E(X)    <= 
### \sum_{i=1}^{k-1}  x_{i+1} (H(x_{i+1}) - H(x_i))
x<-(0:n)/n
Hx<- H(x,alpha=ALPHA,beta=BETA)
ex0<-sum( x[-(n+1)]*( Hx[-1] - Hx[-(n+1)] ) )
ex1<-sum( x[-1]*( Hx[-1] - Hx[-(n+1)] ) )
EH<-(ex0+ex1)/2

### Now find the beta distribution with the same 
### expectation such that pbeta(.05)=.80
rootfunc<-function(x,beta,alpha,EH){
##  x=first parameter, 
## then second parameter is x(1-EH)/EH
	1-beta - pbeta(alpha,x,x*(1-EH)/EH)
}
x<-c((1:higha),1/(1:higha))
px<-pbeta(ALPHA,x,x*(1-EH)/EH)
#plot(x,px,log="x",main=paste("alpha=",ALPHA," beta=",BETA))

if (sign(rootfunc(higha,beta=BETA,alpha=ALPHA,EH=EH))==
	sign(rootfunc(1/higha,beta=BETA,alpha=ALPHA,EH=EH)) 
	){	

mina<- x[px==min(px)]

uroot1<-uniroot(rootfunc,c(1/higha,mina),beta=BETA,
alpha=ALPHA,EH=EH)
a1<-uroot1$root
b1<- a1*(1-EH)/EH 
uroot2<-uniroot(rootfunc,c(mina,higha),beta=BETA,alpha=ALPHA,EH=EH)
a2<-uroot2$root
b2<- a2*(1-EH)/EH 
x<-(0:n)/n
#plot(x,Hx)
#lines(x,pbeta(x,a1,b1),lty=2)
#lines(x,pbeta(x,a2,b2),lty=3)

var1<- var(Hx-pbeta(x,a1,b1))
var2<-var(Hx-pbeta(x,a2,b2))
#print(paste("a1=",a1,"b1=",b1,"a2=",a2,"b2=",b2,
#"var1=",var1,"var2=",var2))
if (var1<var2){ a<-a1 ; b<-b1 }
else { a<-a2 ; b<-b2 }
}
else {
#### else rootfunc(high) and rootfunc(1/higha) 
####have different signs	
uroot<-uniroot(rootfunc,c(1/higha,higha),beta=BETA,
alpha=ALPHA,EH=EH)
a<-uroot$root
b<- a*(1-EH)/EH 
}		
#out<-list(ex0=ex0,ex1=ex1,EH=EH,
#a1=a1,b1=b1,a2=a2,b2=b2,var1=var1,var2=var2,
out<-list(a=a,b=b)
out
}

