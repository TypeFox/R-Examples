prop.comb <-
function(x,n,p=NULL,a=NULL,alternative = c("two.sided", "less", 
    "greater"), conf.level = 0.95)

{




# -------------------------------------------
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


if (!missing(a)) {if (length(a) != length(n)) 
            stop("'a' and 'n' must have the same length")}


#---------------------------------------------





############################
############################
# AUXILIAR FUNCTIONS

test1.MA=function(x,n,p0=0.5,e=5,alternative = c("two.sided", "less", "greater")) {
# Caso de una proporcion y Test
c=1.0/(2.0*n) ; pp=x/n 

dif=pp-p0; part1= dif*dif*n
part2=p0*(1-p0) ;ze0=part1/part2 
 
difc=abs(dif)
if (difc>c) {ze0c=n*(difc-c)*(difc-c)/part2} else {ze0c=0.0} 
 
b1=sqrt((x+0.5)/(n+1.0)); b2=sqrt(p0)
fact=asin(b1)-asin(b2) 
a1=4.0*(n+1.0)*fact*fact 
 
fact1=(x+2.0-(n+4.0)*p0); fact2=(x+2.0)*(n-x+2.0); 
zw2=(fact1*fact1)*(n+4.0)/fact2; 
 
cc1=sqrt(ze0); cc2=sqrt(ze0c); cc3=sqrt(a1); cc4=sqrt(zw2); 
alternative <- match.arg(alternative) 


if(pp<p0){cc1=-cc1;cc2=-cc2;cc3=-cc3;cc4=-cc4} # eran todos positivos

p_ze0=switch(alternative,two.sided =  2.0*pnorm(-abs(cc1)),less = pnorm(cc1),greater = 1-pnorm(cc1))
p_ze0c=switch(alternative,two.sided = 2.0*pnorm(-abs(cc2)),less = pnorm(cc2),greater = 1-pnorm(cc2))
p_a1=switch(alternative,two.sided =   2.0*pnorm(-abs(cc3)),less = pnorm(cc3),greater = 1-pnorm(cc3))
p_zw2=switch(alternative,two.sided =  2.0*pnorm(-abs(cc4)),less = pnorm(cc4),greater = 1-pnorm(cc4))

z=numeric(4)
pvalue=numeric(4)
z[1]=sqrt(ze0); pvalue[1]=p_ze0 # Score (without cc)  
z[2]=sqrt(ze0c); pvalue[2]=p_ze0c # Score (with cc)
z[3]=sqrt(a1); pvalue[3]=p_a1 # Adjusted Arc Sine (a)
z[4]=sqrt(zw2);pvalue[4]=p_zw2 # Adjusted Wald

result=cbind(z,pvalue)
rownames(result)=c("Score (without cc)","Score (with cc)",
"Adjusted Arc Sine", "Adjusted Wald"     )
colnames(result)=c("z value", "p-value")
result} 


ic1.MA=function(x,n,e=5) { 
c=1.0/(2.0*n)
e1=1.0-(e/100.0)/2.0; Z=qnorm(e1); pp=x/n 
 
part1=n/(n+Z*Z); 
part2=(Z*Z)/(4.0*n*n)+pp*(1-pp)/n; 
part3=pp+Z*Z/(2.0*n); 
 
l1ze0=part1*(part3-Z*sqrt(part2)); 
l2ze0=part1*(part3+Z*sqrt(part2)); 
 
part2ccn=(Z*Z)/(4.0*n*n)+(pp-c)*(1-pp+c)/n; 
part2ccp=(Z*Z)/(4.0*n*n)+(pp+c)*(1-pp-c)/n; 
 
l1ze0c=part1*(part3-c-Z*sqrt(part2ccn)); 
l2ze0c=part1*(part3+c+Z*sqrt(part2ccp)); 
 
b1=sqrt((x+0.5)/(n+1.0)); b2=Z/(2.0*sqrt(n+1.0)) 
factl1=asin(b1)-b2; factl2=asin(b1)+b2 
l1a1=sin(factl1)*sin(factl1); l2a1=sin(factl2)*sin(factl2) 
 
fact1=(x+2.0)*(n+2.0-x)/(n+4.0); 
fact2=(x+2.0)-Z*sqrt(fact1); 
fact3=(x+2.0)+Z*sqrt(fact1); 
l1zw2=fact2/(n+4.0); 
l2zw2=fact3/(n+4.0); 
 
part11=(n+Z*Z*Z*Z/53)/(n+Z*Z); 
part21=Z/(n+Z*Z); 
raiz1=n*pp*(1-pp)+Z*Z/4; 
 
l1y0=0.5+part11*(pp-0.5)-part21*sqrt(raiz1); 
l2y0=0.5+part11*(pp-0.5)+part21*sqrt(raiz1); 
 
if (l1ze0<0.0) l1ze0=0.0; if (l2ze0>1.0) l2ze0=1.0; 
if (l1ze0c<0.0) l1ze0c=0.0; if (l2ze0c>1.0) l2ze0c=1.0; 
if (l1a1<0.0) l1a1=0.0; if (l2a1>1.0) l2a1=1.0; 
if (l1zw2<0.0) l1zw2=0.0; if (l2zw2>1.0) l2zw2=1.0; 
if (l1y0<0.0) l1y0=0.0; if (l2y0>1.0) l2y0=1.0; 

result=matrix(nrow=5,ncol=2)
result[1,]=c(l1ze0,l2ze0)
result[2,]=c(l1ze0c,l2ze0c)
result[3,]=c(l1a1,l2a1)
result[4,]=c(l1zw2,l2zw2)
result[5,]=c(l1y0,l2y0)


rownames(result)=c("Score (without cc)","Score (with cc)",
"Adjusted Arc Sine", "Adjusted Wald"    , "Modified Score"  )
colnames(result)=c("lower limit", "upper limit")
result
}

prop.test.MA=function(x,n,p0=0.5,conf.level = 0.95,
alternative = c("two.sided", "less", "greater")) {
e=100*(1-conf.level)
A=test1.MA(x,n,p0,e,alternative)

alternative <- match.arg(alternative) 

B=switch(alternative,two.sided = ic1.MA(x,n,e),less = ic1.MA(x,n,2*e),greater = ic1.MA(x,n,2*e))
B[,1]=switch(alternative,two.sided = B[,1],less = 0,greater=B[,1])
B[,2]=switch(alternative,two.sided = B[,2],less = B[,2],greater=1)

A=rbind(A,c(NA,NA))
result=cbind(B,A)
C=binom.test(x,n,p=p0,alternative=alternative)
C=c(as.numeric(C$conf.int),NA,C$p.value)

result=rbind(result,C)
rownames(result)[6]="Exact binomial test"
list(estimate=x/n,inference=result,alternative=alternative,p0=p0)}

# ESTATISTIC BASED ON ARCSINE TRANSFORMATION
As=function(delta,x2,n2,a,fact1,e1) {
p_vero1=Pvero_dif(delta=delta*sign(a[2]),x=x2,n=n2)
p_vero2=p_vero1 +delta/sign(a[2])
p_vero2[p_vero2<0.0]=0; p_vero2[p_vero2>1.0]=1
fact2=asin(sqrt(p_vero2))-asin(sqrt(p_vero1))
qnorm(e1)^2-4.0*n2[1]*n2[2]*(fact1-fact2)*(fact1-fact2)/(n2[1]+n2[2])}


AM=function(delta,x,n2,a,e1) {
p_vero1=Pvero_dif(delta=delta*sign(a[2]),x=x,n=n2)
p_vero2=p_vero1+delta/sign(a[2])
p_vero2[p_vero2<0.0]=0; p_vero2[p_vero2>1.0]=1
b1=sqrt(p_vero1); b2=sqrt(p_vero2)
b3=sqrt(x[1]/n2[1]); b4=sqrt(x[2]/n2[2])
qnorm(e1)^2-4.0*(n2[1]*(asin(b1)-asin(b3))^2+n2[2]*(asin(b2)-asin(b4))^2)}





ic2.MA=function(x,n,e=5,a=c(-1,1)) { 
e1=1.0-0.5*(e/100.0)
Z=qnorm(e1)^2
aux=FindLambda.MA(x,n,Z=Z,a=a); l1ze0=aux[1]; l2ze0=aux[2]
aux=FindLambda.MA(x,n,Z=Z,a=a,correct=T); l1ze0c=aux[1]; l2ze0c=aux[2]


pp=x/n

I=c(0,0);S=c(0,0)
if(pp[1]==0) {I[1]=1} else {I[1]=0}
if(pp[2]==1) {I[2]=1} else {I[2]=0}
if(pp[1]==1) {S[1]=1} else {S[1]=0}
if(pp[2]==0) {S[2]=1} else {S[2]=0}
 
Z=qnorm(e1)
h=Z^2*(1.0+2.0*I)/4.0
xx=x+h; nn=n+2*h; pp=xx/nn
fact=sum(a^2*(pp*(1-pp))/nn) ; l1zw4=a%*%pp -Z*sqrt(fact) # lower limit

h=Z^2*(1.0+2.0*S)/4.0
xx=x+h; nn=n+2*h; pp=xx/nn
fact=sum(a^2*(pp*(1-pp))/nn); l2zw4=a%*%pp +Z*sqrt(fact)

xx=x+0.5; nn=n+1; b=sqrt(xx/nn)
fact1=asin(b[2])-asin(b[1])

valor=abs(a[2])
aux=uniroot.all(As, c(-1, 1),x=xx,n2=nn,a=a,fact1=fact1,e1=e1)


lim_inf1=min(aux)*valor; lim_sup1=max(aux)*valor

aux=uniroot.all(AM, c(-1, 1),x=xx,n2=nn,a=a,e1=e1)
lim_inf2=min(aux)*valor; lim_sup2=max(aux)*valor

xx=x+Z^2/4; nn=n+Z^2/2

aux=FindLambda.MA(x=xx,n=nn,Z=Z^2,correct=F,a);
lim_inf3=min(aux)*valor; lim_sup3=max(aux)*valor

aux=FindLambda.MA(x=xx,n=nn,Z=Z^2,correct=T,a);
lim_inf4=min(aux)*valor; lim_sup4=max(aux)*valor

l1ze3=lim_inf3;l2ze3=lim_sup3;
l1ze3c=lim_inf4;l2ze3c=lim_sup4;

l1ze3=max(-1.0,l1ze3);l2ze3=min(1.0,l2ze3)
l1ze3c=max(-1.0,l1ze3c);l2ze3c=min(1.0,l2ze3c)

l1ae1=lim_inf1;l2ae1=lim_sup1;
l1aem1=lim_inf2;l2aem1=lim_sup2;

l1ze0=max(-1.0,l1ze0);l2ze0=min(1.0,l2ze0)
l1ze0c=max(-1.0,l1ze0c);l2ze0c=min(1.0,l2ze0c)
l1ae1=max(-1.0,l1ae1);l2ae1=min(1.0,l2ae1)
l1aem1=max(-1.0,l1aem1);l2aem1=min(1.0,l2aem1)
l1zw4=max(-1.0,l1zw4);l2zw4=min(1.0,l2zw4)

result=cbind(c(l1ze0,l1ze0c,l1ze3,l1ze3c,l1zw4,l1ae1,l1aem1),
c(l2ze0,l2ze0c,l2ze3,l2ze3c,l2zw4,l2ae1,l2aem1))

rownames(result)=c("Score (without cc)","Score (with cc)","Adjusted Score (without cc)", 
"Adjusted Score (with cc)","Adjusted Wald","Adjusted Arc Sine","Adjusted-M Arc Sine"  )
colnames(result)=c("lower limit", "upper limit")
result}



test2.MA=function(x,n,a,lambda=0,e=5,
alternative = c("two.sided", "less", "greater")) {

alternative <- match.arg(alternative) 
ze0=FindZ.MA(x,n,correct=FALSE,a=a,lambda=lambda,plot=TRUE)
ze0c=FindZ.MA(x,n,correct=TRUE,a=a,lambda=lambda,plot=TRUE)

e1=1.0-(e/100.0)/2.0;Z=qnorm(e1)

xx=x+Z^2/4;nn=n+Z^2/2
p=xx/nn; d=a%*%p
#pp=sum(xx)/sum(nn);qq=1.0-pp;
pp=Pvero_dif(delta=lambda/a[2],x=xx,n=nn); qq=1-pp

fact1=a[1]^2*pp*qq/nn[1]+a[2]^2*(pp+lambda)*(1-pp-lambda)/nn[2]; fact2=nn[1]*nn[2]+nn[1]+nn[2];
if (nn[1]==nn[2]) {c=2.0/fact2} else {c=1.0/fact2}

ze3=sqrt((d-lambda)^2/fact1) #para cualquier valor de lambda (no especificamente 0)

if (d-lambda>0) {dc=d-lambda} else {dc=-d+lambda}
if (dc>c) {e3c=(dc-c)^2/fact1} else {e3c=0.0}
ze3c=sqrt(e3c)

ze0=sqrt(ze0);ze0c=sqrt(ze0c)

pp=x/n

I=c(0,0);S=c(0,0)
if(pp[1]==0) {I[1]=1} else {I[1]=0}
if(pp[2]==1) {I[2]=1} else {I[2]=0}
if(pp[1]==1) {S[1]=1} else {S[1]=0}
if(pp[2]==0) {S[2]=1} else {S[2]=0}

dif=a%*%pp
 
if(dif>lambda) 
{h=Z^2*(1.0+2.0*I)/4.0} else
{h=Z^2*(1.0+2.0*S)/4.0}

xx=x+h; nn=n+2*h; pp=xx/nn
fact=sum(a^2*(pp*(1-pp))/nn) 
zw4=(a%*%pp -lambda)/sqrt(fact) 


xx=x+0.5; nn=n+1; b=sqrt(xx/nn)
fact1=asin(b[2])-asin(b[1])

aux=uniroot.all(AM, c(-1, 1),x=xx,n2=nn,a=a,e1=e1)
p_vero1=Pvero_dif(delta=lambda/a[2],x=xx,n=nn)
p_vero2=p_vero1+lambda/a[2]
if(p_vero2<0.0) p_vero2=0.0; if(p_vero2>1.0) p_vero2=1.0

b=c(b,sqrt(p_vero1),sqrt(p_vero2))
fact2=asin(b[2])-asin(b[1])-asin(b[4])+asin(b[3]);
ae1=4.0*nn[1]*nn[2]*fact2^2/(nn[1]+nn[2])

ae1m=4.0*(nn[1]* ( asin(b[1]) - asin(b[3]) )^2 + nn[2]* (asin(b[2]) - asin(b[4]) ) ^2)

cc3=sqrt(ae1);cc5=sqrt(ae1m);cc4=zw4;


zs=c(ze0,ze0c,ze3,ze3c,cc4,cc3,cc5);zs=abs(zs)
pvalue=numeric(7)

if(a%*%(x/n)<lambda){zs=(-zs)} # eran todos positivos


pvalue=switch(alternative,two.sided =  2.0*pnorm(-abs(zs)),less = pnorm(zs),greater = 1-pnorm(zs))

result=cbind(zs,pvalue)
rownames(result)=c("Score (without cc)","Score (with cc)","Adjusted Score (without cc)", 
"Adjusted Score (with cc)","Adjusted Wald","Adjusted Arc Sine", "Adjusted-M Arc Sine")
colnames(result)=c("z value", "p-value")
result} 



prop.test2.MA=function(x,n,lambda,conf.level=0.95,
alternative = c("two.sided", "less", "greater")) {
a=c(-1,1)
e=100*(1-conf.level)

A=test2.MA(x,n,a,lambda,e,alternative)

alternative <- match.arg(alternative) 
B=switch(alternative,two.sided = ic2.MA(x,n,e,a),less = ic2.MA(x,n,e=2*e,a),greater = ic2.MA(x,n,e=2*e,a))
B[,1]=switch(alternative,two.sided = B[,1],less = -1,greater=B[,1])
B[,2]=switch(alternative,two.sided = B[,2],less = B[,2],greater=1)

result=cbind(B,A)

list(estimate=x/n,a=a,difference=as.numeric(a%*%(x/n)),inference=result,alternative=alternative,lambda=lambda)
}



testComb2.MA=function(x,n,a,lambda=0,e=5,
alternative = c("two.sided", "less", "greater")) {
alternative = match.arg(alternative) 

prop=x/n; L=a%*%prop

ze0=FindZ.MA(x,n,correct=FALSE,a=a,lambda=lambda)
ze0c=FindZ.MA(x,n,correct=TRUE,a=a,lambda=lambda)

xx=x+1;nn=n+2; p=xx/nn; L=a%*%p; d=L-lambda

fact=sum(a^2*p*(1-p)/nn); zw2=d^2/fact;

e1=1.0-0.5*(e/100);Z=qnorm(e1)

if (a[1]>0) {s1=1.0} else {s1=-1.0}
if (a[2]>0) {s2=1.0} else {s2=-1.0}

pp11=0.5*(1.0+s1);pp12=0.5*(1.0+s2)
pp21=0.5*(1.0-s1);pp22=0.5*(1.0-s2)

if (pp11==prop[1]) {I1=1} else {I1=0}
if (pp12==prop[2]) {I2=1} else {I2=0}
if (pp21==prop[1]) {S1=1} else {S1=0}
if (pp22==prop[2]) {S2=1} else {S2=0}

if (a%*%prop>lambda)
{h=c(1+I1*2,1+I2*2)} else {h=c(1+S1*2,1+S2*2)}

xxx=x+Z^2*h/4; nnn=n+Z^2*h/2; pd=xxx/nnn; Ld=a%*%pd; dd=Ld-lambda
factd=sum(a^2*pd*(1-pd)/nnn)
zw4=dd^2/factd


# caso de LR1*/
xm=x+0.5; nm=n+1;pm=xm/nm

#Calcular la p_verosimilitud con los datos incrementados*/

ad1=sum(xm)
alpha=lambda/a[2];bet=-a[1]/a[2]
c0=alpha*(1.0-alpha)*xm[1]
c1=bet*ad1-nm[1]*alpha*(1.0-alpha)-alpha*bet*(nm[2]+2.0*xm[1])
c2=bet*(alpha*(nm[2]+2*nm[1])-(nm[1]+xm[2])-bet*(nm[2]+xm[1]))
c3=bet*bet*(nm[1]+nm[2])

A=-4.5*c3*(c1*c2-3.0*c0*c3)+c2*c2*c2  #        /*Recordar que -A*/
B=c2*c2-3.0*c1*c3

if (A>=0)  {u=sqrt(B)} else {u=-sqrt(B)}
w1= A/(u^3); if (w1>1.0) w1=1.0; if (w1<(-1.0)) w1=-1.0
	
fhi=(pi+ acos(w1))/3
pvm1=(-1*c2+2*cos(fhi)*u)/(3*c3) ;if (pvm1>1) pvm1=1; if (pvm1<0) pvm1=0
pvm2=alpha+bet*pvm1;if (pvm2>1) pvm2=1;if (pvm2<0) pvm2=0; 



gk1=log(pvm1/pm[1]);gk2=log((1-pvm1)/(1-pm[1]));gk3=log(pvm2/pm[2]);gk4=log((1-pvm2)/(1-pm[2]));
lr1=-2*(xm[1]*gk1+(nm[1]-xm[1])*gk2+xm[2]*gk3+(nm[2]-xm[2])*gk4)

cc=c(sqrt(ze0),sqrt(ze0c),sqrt(zw2),sqrt(zw4),sqrt(lr1))
if(a%*%(x/n)<lambda){cc=-cc} # eran todos positivos


pvalor=switch(alternative,two.sided =  2.0*pnorm(-abs(cc)),less = pnorm(cc),greater = 1-pnorm(cc))

result=cbind(cc,pvalor)
rownames(result)=c("Score (without cc)","Score (with cc)",
"Adjusted Wald(a)", "Adjusted Wald(b)"  ,"Adjusted likelihood ratio"  )
colnames(result)=c("z value", "p-value")
result
}





icComb2.MA=function(x,n,e=5,a) { 
prop=x/n; e1=1-(e/100)/2; Z=qnorm(e1)^2
limites=c(sum(a[a<0]),sum(a[a>0]))
lze0=FindLambda.MA(x,n,Z=Z,correct=FALSE,a=a)
lze0c=FindLambda.MA(x,n,Z=Z,correct=TRUE,a=a)

xx=x+1;nn=n+2; p=xx/nn; L=a%*%p; fact=sum(a^2*p*(1-p)/nn)
Z=qnorm(e1)

lzw2=c(L-Z*sqrt(fact),L+Z*sqrt(fact))
if (a[1]>0) {s1=1.0} else {s1=-1.0}
if (a[2]>0) {s2=1.0} else {s2=-1.0}
pp11=0.5*(1.0+s1);pp12=0.5*(1.0+s2)
pp21=0.5*(1.0-s1);pp22=0.5*(1.0-s2)

if (pp11==prop[1]) {I1=1} else {I1=0}
if (pp12==prop[2]) {I2=1} else {I2=0}
if (pp21==prop[1]) {S1=1} else {S1=0}
if (pp22==prop[2]) {S2=1} else {S2=0}

#limite inferior
h=c(1+I1*2,1+I2*2)
xxx=x+Z^2*h/4; nnn=n+Z^2*h/2; pd=xxx/nnn; Ldi=a%*%pd
factdi=sum(a^2*pd*(1-pd)/nnn)

#limite superior
h=c(1+S1*2,1+S2*2)
xxx=x+Z^2*h/4; nnn=n+Z^2*h/2; pd=xxx/nnn; Lds=a%*%pd
factds=sum(a^2*pd*(1-pd)/nnn)

lzw4=c(Ldi-Z*sqrt(factdi),Lds+Z*sqrt(factds))

# /*caso de LR1*/
xm=x+0.5;nm=n+1;pm=xm/nm
llr=uniroot.all(RazonVero, limites,x=xm,n2=nm,a=a,Z=Z)


result=rbind(lze0,lze0c,lzw2,lzw4,llr)
rownames(result)=c("Score (without cc)","Score (with cc)","Adjusted Wald (a)",
"Adjusted Wald (b)", "Adjusted likelihood ratio"  )
colnames(result)=c("lower limit", "upper limit")

result[result[,1]<limites[1]]=limites[1]; result[result[,1]>limites[2]]=limites[2]
result}



prop.Comb2.MA=function(x,n,a,lambda,conf.level=0.95,
alternative = c("two.sided", "less", "greater")) {

e=100*(1-conf.level)
A=testComb2.MA(x,n,a,lambda,e,alternative)

alternative =match.arg(alternative) 

B=switch(alternative,two.sided = icComb2.MA(x,n,e,a),less = icComb2.MA(x,n,e=2*e,a),
greater = icComb2.MA(x,n,e=2*e,a))

B[,1]=switch(alternative,two.sided = B[,1],less = sum(a[a<0]),greater=B[,1])
B[,2]=switch(alternative,two.sided = B[,2],less = B[,2],greater=sum(a[a>0]))

result=cbind(B,A)

list(estimate=x/n,a=a,L=as.numeric(a%*%(x/n)),inference=result,alternative=alternative,lambda=lambda)
}

testCombi.MA=function(x,n,a,lambda=0,e=5,
alternative = c("two.sided", "less", "greater")) {

alternative = match.arg(alternative) 
ze0=FindZ.MA(x,n,correct=FALSE,a=a,lambda=lambda,plot=TRUE)
ze0c=FindZ.MA(x,n,correct=TRUE,a=a,lambda=lambda,plot=TRUE)

e1=1.0-0.5*(e/100); Z=qnorm(e1) 
fact1=0.5*Z^2; pr=x/n ;L=a%*%pr

K=length(pr); prop_sa=numeric(K)
s=rep(1,K);s[a<0]=-1 
prop_sa=0.5*(1+s); prop_sb=0.5*(1-s)

h=numeric(K)
if (L>lambda){ii=prop_sa==pr; h[ii]=fact1*(1+K)/K; h[!ii]=fact1/K}
else {ii=prop_sb==pr; h[ii]=fact1*(1+K)/K; h[!ii]=fact1/K}

xx=x+h;nn=n+2*h;pp=xx/nn; L=a%*%pp
denom=sum(a^2*pp*(1-pp)/nn); denom2=sum(a^2/n)

zw3=sqrt((L-lambda)^2/denom)

B=sum(a); denom3=denom2-(B-2*lambda)^2/sum(n)
pa0=4.0*(a%*%pr-lambda)^2/denom3; zpa0=sqrt(pa0)

z=numeric(4)
z[1]=sqrt(ze0); z[2]=sqrt(ze0c);z[3]=zw3;z[4]=zpa0

if(a%*%(x/n)<lambda){z=-z} # eran todos positivos

pvalue=switch(alternative,two.sided = 2.0*pnorm(-abs(z)),less = pnorm(z),greater = 1-pnorm(z))

result=cbind(z,pvalue)
rownames(result)=c("Score (without cc)","Score (with cc)",
"Adjusted Wald", "Peskun "); colnames(result)=c("z value", "p-value")
result}





icCombi.MA=function(x,n,e=5,a) { 
e1=1.0-0.5*(e/100.0);Z=qnorm(e1)^2

aux=FindLambda.MA(x,n,Z=Z,a=a); l1ze0=aux[1]; l2ze0=aux[2]
aux=FindLambda.MA(x,n,Z=Z,correct=T,a=a); l1ze0c=aux[1]; l2ze0c=aux[2]

fact1=0.5*Z; pr=x/n ;L=a%*%pr
K=length(pr); prop_sa=numeric(K)
s=rep(1,K);s[a<0]=-1 
prop_sa=0.5*(1+s); prop_sb=0.5*(1-s)

h1=numeric(K);h2=numeric(K)
ii=prop_sa==pr; h1[ii]=fact1*(1+K)/K; h1[!ii]=fact1/K
ii=prop_sb==pr; h2[ii]=fact1*(1+K)/K; h2[!ii]=fact1/K

xx1=x+h1;nn1=n+2*h1;pp1=xx1/nn1; L1=a%*%pp1
xx2=x+h2;nn2=n+2*h2;pp2=xx2/nn2; L2=a%*%pp2
denom1=sum(a^2*pp1*(1-pp1)/nn1)
denom2=sum(a^2*pp2*(1-pp2)/nn2)

lim=matrix(ncol=2,nrow=4)
lim[1,]=c(l1ze0,l2ze0); lim[2,]=c(l1ze0c,l2ze0c)
lim[3,]=c(L1-qnorm(e1)*sqrt(denom1),L2+qnorm(e1)*sqrt(denom2))

B=sum(a);sumn=sum(n); c1=sumn/(sumn+Z);
L2=a%*%pr; c2=L2+B*Z/(2*sumn);c3=(B-2*L2)^2/sumn;
c4=sum(a^2/n)*(sumn+Z)/sumn;c5=c4-c3;

lim[4,]=c(c1*(c2-qnorm(e1)*sqrt(c5)*0.5),c1*(c2+qnorm(e1)*sqrt(c5)*0.5))

B_neg=sum(a[a<0]);B_pos=sum(a[a>0])

ii=lim[,1]<B_neg; lim[ii,1]=B_neg
ii=lim[,2]>B_pos; lim[ii,2]=B_pos

rownames(lim)=c("Score (without cc)","Score (with cc)","Adjusted Wald",
"Peskun ")
colnames(lim)=c("lower limit", "upper limit")
lim}



prop.Combi.MA=function(x,n,a,lambda,conf.level=0.95,
alternative = c("two.sided", "less", "greater")) {

e=100*(1-conf.level)

A=testCombi.MA(x,n,a,lambda,e,alternative)

alternative =match.arg(alternative) 

B=switch(alternative,two.sided = icCombi.MA(x,n,e,a),less = icCombi.MA(x,n,e=2*e,a),
greater = icCombi.MA(x,n,e=2*e,a))

B[,1]=switch(alternative,two.sided = B[,1],less = sum(a[a<0]),greater=B[,1])
B[,2]=switch(alternative,two.sided = B[,2],less = B[,2],greater=sum(a[a>0]))

result=cbind(B,A)

list(estimate=x/n,a=a,L=as.numeric(a%*%(x/n)),inference=result,alternative=alternative,
lambda=lambda)
}

Pvero_dif=function(delta,x,n){
	a1=sum(x); a=sum(n)
	b=-(a+a1-(a+n[1])*delta)
	c=-delta*(-n[1]*delta+2.0*x[1]+sum(n))+a1
	d=x[1]*delta*(1-delta)

	v1=b/(3.0*a);v2=v1^3;v3=(b*c)/(6.0*a^2);v4= d/(2.0*a)
	v=v2-v3+v4
      u= sqrt( b^2/(9.0*a^2)-c/(3*a))

      ii=v<0; u[ii]=-u[ii]     
	w1= v/(u^3)
       w1[w1>1]=1.0; w1[w1<(-1.0)]=-1.0
	w2=acos(w1); w=(acos(-1.0)+w2)/3.0

	pvero=2.0*u*cos(w)-b/(3*a)
	pvero[pvero>1.0]=1.0; pvero[pvero<0.0]=0 
	pvero}





###########################
###########################

# END OF AUXILIAR FUNCTIONS

############################
############################


OK=complete.cases(x, n); x=x[OK]; n=n[OK]
k=length(x)

recomen=""

##############################################################

#               A PROPORTION


if (k==1) {
if (missing(p)) p=0.5
else if (p<=0 | p>=1) stop("value of 'p' must be in (0,1)")


resultado=prop.test.MA(x,n,p0=p,conf.level,alternative)

# Recomendation ____

if (conf.level==0.90 & sum(n)>=100) {recomen="Adjusted Arc Sine"}
else if (conf.level!=0.90 & sum(n)<=80) {recomen="Score (with cc)"}
else {recomen= "Adjusted Wald"} }




##############################################################

#               DIFFERENCE OF PROPORTIONS

else if(k==2 & missing(a)) {

  if (missing(p)) p=0
  else if (p<=-1 | p>=1) stop("value of 'p' must be in (-1,1)")


  resultado=prop.test2.MA(x=x,n=n,lambda=p,conf.level,alternative)
  a=resultado$a
  
# Recommendation ____

if (p==0) {
if (conf.level==0.90) recomen="Adjusted Score (without cc)"
if (conf.level!=0.90 & n[1]==n[2]) recomen="Adjusted Score (with cc)"
if (conf.level!=0.90 & n[1]!=n[2]) recomen="Score (with cc)"}
else {
if (conf.level==0.99) recomen="Adjusted-M Arc Sine"
if (conf.level!=0.99 & sum(n)>=200 & n[1]==n[2]) recomen="Adjusted Wald"
if (conf.level!=0.99 & (sum(n)<200 | n[1]!=n[2])) recomen="Adjusted Arc Sine"
}}



##############################################################

#               COMBINATION OF TWO PROPORTIONS

else if(k==2 & !missing(a)) {
  if (missing(p)) p=0.5*(sum(a[a<0])+ sum(a[a>0]) )
  else {
   ai=sum(a[a<0]); as=sum(a[a>0])
   if (p<=ai | p>=as) {
   cat(paste("a_neg=",ai,": sum of negative coefficients",sep=""),"\n")
cat(paste("a_pos=",as,": sum of positive coefficients",sep=""),"\n")

stop("value of 'p' must be in (a_neg,a_pos)")
}}


  resultado=prop.Comb2.MA(x=x,n=n,a=a,lambda=p,conf.level,alternative)
 
# Recomendation ____

if (conf.level!=0.90 &  n[1]<60 & n[2]<60) recomen="Adjusted Wald(a)"
if (conf.level!=0.99) recomen="Adjusted Wald(b)"
if (conf.level==0.99) recomen="Adjusted Likelihood Ratio"

 }



##############################################################

#               COMBINATION OF MORE THAN TWO PROPORTIONS


else if(k>2) {

if (missing(a)) stop ("value of 'a' is lost")

if (missing(p)) p=0.5*(sum(a[a<0])+ sum(a[a>0]) )
else {
   ai=sum(a[a<0]); as=sum(a[a>0])
  
  if (p<=ai | p>=as) {


cat(paste("a_neg=",ai,": sum of negative coefficients",sep=""),"\n")
cat(paste("a_pos=",as,": sum of positive coefficients",sep=""),"\n")
stop("value of 'p' must be in (a_neg,a_pos)")}}


resultado=prop.Combi.MA(x=x,n=n,a=a,lambda=p,conf.level,alternative)

if (k==3 & max(n)>10) recomen="Score (with cc)"
else if (k>3 & max(n)>10) recomen="Score (without cc)"
else if ( max(n)<=10) recomen="Adjusted Wald"

}


if (conf.level!=0.90 & conf.level!=0.95 &conf.level!=0.99) {recomen=""}

if (recomen=="") recomen="No recommendation"


resultado$k=k; resultado$x=x; resultado$n=n;
resultado$a=a; resultado$p=p; resultado$conf.level=conf.level;
resultado$recomendation=recomen
object=resultado
class(object) ="prop.comb"
object
}
