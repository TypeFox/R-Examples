prop.RR <-
function(x,n,rho=NULL,alternative = c("two.sided", "less", 
    "greater"), conf.level = 0.95)
{

# -----------------------------------------

if (length(x) != length(n)) 
            stop("'x' and 'n' must have the same length")
   
 if (any(n <= 0)) 
  stop("elements of 'n' must be positive")
 if (any(x < 0)) 
   stop("elements of 'x' must be nonnegative")
 if (any(x > n)) 
 stop("elements of 'x' must not be greater than those of 'n'")


 xr <- round(x)
 if (any(is.na(x) | (x < 0)) || max(abs(x - xr)) > 1e-07) 
 stop("'x' must be nonnegative and integer")
 x=xr

 nr <- round(n)
 if (any(is.na(n) | (n < 1)) || max(abs(n - nr)) > 1e-07) 
 stop("'n' must be a positive integer")
 n=nr

if (!missing(rho)) {
if (rho<=0) 
stop("'rho' must be positive")
}
# ------------------------

##############################
###   AUXILIAR FUNCTIONS


ic.RR.MA=function(x,n,e=5) { 
e1=1-0.5*(e/100); Z=qnorm(e1)
xx=x+0.5;nn=n+1;p=xx/nn
xx1=xx[1];xx2=xx[2];nn1=nn[1];nn2=nn[2];p1=p[1];p2=p[2]

c1=sum(nn)*xx[1]*xx[2]+Z^2*(nn[1]*xx[1]+nn[2]*xx[2]-2*xx[1]*xx[2])*0.5
c2=sum(nn)^2*xx[1]*xx[2]*(sum(xx)-sum(nn)*p[1]*p[2])
c3=Z*(nn[2]*xx[2]-nn[1]*xx[1])*0.5
c4=xx[1]*(sum(nn)*nn[2]*p[1]-Z^2*(nn[2]-xx[1]))

aux=Z*sqrt(c2+c3^2); zab1=c(c1-aux,c1+aux);zab1=zab1/c4

casoesp=sum(nn)*nn[2]*xx[1]/(nn[1]*(nn[2]-xx[1]))

lzab1=numeric(2)

if(Z^2==casoesp){
   s1=(nn2-xx1)*p2/p1+(nn1-xx2)
   s2=(nn1+nn2+nn1)*(nn2-xx1)*p2/p1+nn1*(nn1-xx2)
   lzab1=c((1.0-(nn1+nn2)*s1/s2)*p2/p1,Inf)}

else{ 
	if ((zab1[1]>=xx2/(nn1+nn2-xx1))&(zab1[1]<p2/p1)) {lzab1[1]=zab1[1]}
  	else {if ((zab1[2]>=xx2/(nn1+nn2-xx1) ) & (zab1[2]<p2/p1)) {lzab1[1]=zab1[2]} 
   			else {mi1=1/(nn2*p1^2+Z^2);mi2=xx2*p1+0.5*Z^2;mi3=0.25*Z^2+xx2*(p1-p2);
	 		lzab1[1]=mi1*(mi2-Z*sqrt(mi3))}}
	
     if ( (zab1[2]<=(nn1+nn2-xx2)/xx1 ) & (zab1[2]>p2/p1)) {lzab1[2]=zab1[2]}
     else {if ( (zab1[1]<=(nn1+nn2-xx2)/xx1 )&(zab1[1]>p2/p1) ) {lzab1[2]=zab1[1]} 
			else {ms1=1/(nn1*p1^2);ms2=xx1*p2+0.5*Z^2;ms3=0.25*Z^2+xx1*(p2-p1);
			lzab1[2]=ms1*(ms2+Z*sqrt(ms3))}}
    }

prop=x/n

if(prop[1]==0) {I1=1} else {I1=0}; if(prop[2]==1) {I2=1} else {I2=0}
if(prop[1]==1) {S1=1} else {S1=0}; if(prop[2]==0) {S2=1} else {S2=0}

h11=1+I1*2;h12=1+I2*2; h21=1+S1*2;h22=1+S2*2

xxi1=x[1]+Z*Z*h11/4;xxi2=x[2]+Z*Z*h12/4
nni1=n[1]+Z*Z*h11/2;nni2=n[2]+Z*Z*h12/2
p1di=xxi1/nni1; p2di=xxi2/nni2

xxs1=x[1]+Z^2*h21/4;xxs2=x[2]+Z^2*h22/4
nns1=n[1]+Z^2*h21/2;nns2=n[2]+Z^2*h22/2
p1ds=xxs1/nns1;p2ds=xxs2/nns2

casoesp2i=nni1*p1di/(1.0-p1di)

if (casoesp2i==Z^2) {l1zw4=0.5*( p2di/p1di-nni1*(1.0-p2di)/( nni2*(1.0-p1di) ) )}
else{
fact1i=(nni1-xxi1)/(nni1*xxi1)+(nni2-xxi2)/(nni2*xxi2)-Z^2*(nni1-xxi1)*(nni2-xxi2)/(nni1*xxi1*nni2*xxi2);
fact2i=(p2di/p1di)/(1-Z^2*(nni1-xxi1)/(nni1*xxi1) )
l1zw4=fact2i*(1-Z*sqrt(fact1i))}

casoesp2s=nns1*p1ds/(1-p1ds)

if(casoesp2s==Z^2) {l2zw4=Inf} else{
fact1s=(nns1-xxs1)/(nns1*xxs1)+(nns2-xxs2)/(nns2*xxs2)-Z*Z*(nns1-xxs1)*(nns2-xxs2)/(nns1*xxs1*nns2*xxs2)
fact2s=p2ds/(p1ds*(1-Z^2*(nns1-xxs1)/(nns1*xxs1)))
l2zw4=fact2s*(1+Z*sqrt(fact1s))}
lzw4=c(l1zw4,l2zw4)

R_m=p2/p1;fact=(nn1-xx1)/(nn1*xx1)+(nn2-xx2)/(nn2*xx2)
l1lw1=R_m*exp(-Z*sqrt(fact)); l2lw1=R_m*exp(+Z*sqrt(fact))

llw1=c(l1lw1,l2lw1)

#######################################################

# limite donde buscar la solucion acotada por min de limites inferiores y 
# max de limites superiores

ls=max(l2lw1,l2zw4,lzab1[2]); li=min(l1lw1,l1zw4,lzab1[1])


Score_RR_MA=function (rho,x,n2,e) {
e1=1-0.5*(e/100)
p=x/n2; p1=Pvero_ratio(x,n2,rho); p2=rho*p1
qnorm(e1)^2-(p[2]-rho*p[1])^2/(rho^2*p1*(1-p1)/n2[1]+p2*(1-p2)/n2[2])
}



aux=uniroot.all(Score_RR_MA, c(li,ls),x=x,n2=n,e=e)
lze0=c(min(aux),max(aux))


AM2=function(rho,x,n2,e) {
e1=1-0.5*(e/100)
p_vero1=Pvero_ratio(x,n2,rho); p_vero2=p_vero1*rho
fact2=asin(sqrt(p_vero2))-asin(sqrt(x[2]/n2[2]))
fact1=asin(sqrt(p_vero1))-asin(sqrt(x[1]/n2[1]))
qnorm(e1)^2-4*(n2[1]*fact1^2+ n2[2]*fact2^2)}


aux=uniroot.all(AM2,c(li,ls),x=x+0.5,n2=n+1,e=e)
lal1m=c(min(aux),max(aux))

aux=uniroot.all(AM2,c(li,ls),x=x+1,n2=n+2,e=e)
lal2m=c(min(aux),max(aux))

result=rbind(  lze0, lzab1 ,lzw4,llw1,lal1m,lal2m)
i1=result[,1]<0;result[i1,1]=0


rownames(result)=c(
"Score (without cc)",
"Adjusted Score approx.",
"Adjusted Wald",
"Adjusted log transfor.", 
"Adjusted-M Arc Sine (a)",
"Adjusted-M Arc Sine (b)"   )

colnames(result)=c("lower limit", "upper limit")
result}


test.RR.MA=function(x,n,rho,e=5,
alternative = c("two.sided", "less", "greater")) {
alternative <- match.arg(alternative) 


ze0=FindZ.MA(x=x,n=n,a=c(-rho,1))
e1=1-(e/100)/2;Z=qnorm(e1)

xx1=x[1]+0.5;xx2=x[2]+0.5;nn1=n[1]+1.0;nn2=n[2]+1.0;

p1=xx1/nn1;p2=xx2/nn2;

if ((rho>=xx2/(nn1+nn2-xx1))&(rho<=(nn1+nn2-xx2)/xx1)){
m1=(nn1+nn2)*nn1*nn2;m2=xx2+rho*xx1;
m3=nn1-xx2+rho*(nn2-xx1);zab1=m1*(p2-rho*p1)*(p2-rho*p1)/(m2*m3)}
else {
if (rho<xx2/(nn1+nn2-xx1)) {zab1=nn2*(p2-rho*p1)^2/(rho*(1-rho))}
	  else  {zab1=nn1*(p2-rho*p1)^2/(rho-1.0)}}

prop=x/n

if(prop[1]==0) {I1=1} else {I1=0}; if(prop[2]==1) {I2=1} else {I2=0}
if(prop[1]==1) {S1=1} else {S1=0}; if(prop[2]==0) {S2=1} else {S2=0}

coc=prop[2]/prop[1]
if(coc>rho){h=c(1+I1*2,1+I2*2)} else {h=c(1+S1*2,1+S2*2)}



xx=x+Z^2*h/4; nn=n+Z^2*h/2; pd=xx/nn
fact1=pd[2]-rho*pd[1];fact2=rho^2*pd[1]*(1-pd[1])/nn[1]+ pd[2]*(1-pd[2])/nn[2]  
zw4=fact1^2/fact2

R_m=p2/p1;fact=1/xx1+1/xx2-(nn1+nn2)/(nn1*nn2);lw1=log(R_m/rho)^2/fact


AM3=function(rho,x,n2) {
p_vero1=Pvero_ratio(x,n2,rho); p_vero2=p_vero1*rho
fact2=asin(sqrt(p_vero2))-asin(sqrt(x[2]/n2[2]))
fact1=asin(sqrt(p_vero1))-asin(sqrt(x[1]/n2[1]))
4*(n2[1]*fact1^2+ n2[2]*fact2^2)}

alm=c(AM3(rho,x=x+0.5,n2=n+1),AM3(rho,x=x+1,n2=n+2))

z=sqrt(c(ze0,zab1,zw4,lw1,alm))

p=x/n; if(p[2]<rho*p[1]){z=-z} # eran todos positivos
pvalor=switch(alternative,two.sided = 2*pnorm(-abs(z)),less = pnorm(z),greater=1-pnorm(z))

result=cbind(z,pvalor)
rownames(result)=c("Score (without cc)","Adjusted score approx.",
"Adjusted Wald", "Adjusted log transfor." ,"Adjusted-M arc sine (a)","Adjusted-M arc sine (b)"  )
colnames(result)=c("z value", "p-value")
result
}


RR.MA=function(x,n,rho,conf.level=0.95,
alternative = c("two.sided", "less", "greater")) {

e=100*(1-conf.level)

A=test.RR.MA(x,n,rho,e,alternative)
alternative =match.arg(alternative) 

B=switch(alternative,two.sided = ic.RR.MA(x,n,e),less = ic.RR.MA(x,n,2*e),
greater = ic.RR.MA(x,n,2*e))

B[,1]=switch(alternative,two.sided = B[,1],less = 0 ,greater=B[,1])
B[,2]=switch(alternative,two.sided = B[,2],less = B[,2],greater=Inf)

result=cbind(B,A)

p=x/n
list(estimate=p,RR=p[2]/p[1],inference=result,alternative=alternative,rho=rho)
}




##############################
###   END OF AUXILIAR FUNCTIONS
###############################




recomen=""
OK=complete.cases(x, n); x=x[OK]; n=n[OK]
k=length(x)

if (missing(rho)) rho=1

resultado=RR.MA(x,n,rho,conf.level,alternative) 

resultado$x=x; resultado$n=n;
resultado$rho=rho; resultado$conf.level=conf.level


if (conf.level==0.99) recomen="Adjusted-M Arc Sine (b)"
else if (conf.level!=0.99 & (rho<0.1 | rho>10) ) recomen="Adjusted-M Arc Sine (b)"
else if (conf.level!=0.99 & (rho>=0.1 & rho<=10) ) recomen="Adjusted Score approx."
else if (conf.level==0.95 & sum(n)>200 & (rho>=0.1 & rho<=10) ) recomen="Adjusted log transfor."

else if ((rho>=0.1 & rho<=10)) recomen="Adjusted Wald"
else if (sum(n)>200 & conf.level==0.95 & (rho>=0.2 & rho<=5))  recomen=="Score (without cc)" 

if (conf.level!=0.90 & conf.level!=0.95 &conf.level!=0.99) {recomen=""}
if (recomen=="") recomen="No recommendation"


resultado$recomendation=recomen
object=resultado
class(object) ="prop.RR"
object
}
