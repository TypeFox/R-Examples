#File:		QM.R
#Version:	1.1-2
#Author:	Juha Karvanen
#Date:		11 October 2005

data2normpoly4<-function(data)
{
#data is vector
lmoms<-Lmoments(data)
params<-lmom2normpoly4(lmoms)
params
}

data2normpoly6<-function(data)
{
#data is vector
lmoms<-Lmoments(data)
params<-lmom2normpoly6(lmoms)
params
}

lmom2normpoly4<-function(lmom)
{
	L1<-lmom[1];
	L2<-lmom[2];
	L3<-lmom[3];
	L4<-lmom[4];

	a1<-L1-3*L2+5*L3+24.46947735493778*L4;
	a2<-6*L2-30*L3-48.93895470987556*L4;
	a3<-30*L3;
	a4<-14.457006455887834*L4;

	param<-c(a1,a2,a3,a4);
}

lmom2normpoly6<-function(lmom)
{
	L1<-lmom[1];
	L2<-lmom[2];
	L3<-lmom[3];
	L4<-lmom[4];
	L5<-lmom[5];
	L6<-lmom[6];

	a1<-L1-3*L2+5*L3-7*L4+9*L5+88.36715690214288*L6; 
	a2<-6*L2-30*L3+84*L4-180*L5-373.2962367553012*L6; 
	a3<-30*L3-210*L4+810*L5+589.6857688530462*L6; 
	a4<-140*L4-1260*L5-393.1238459020308*L6; 
	a5<-630*L5; 
	a6<-40.595671272636515*L6;

	param<-c(a1,a2,a3,a4,a5,a6);
}




qnormpoly<-function(cp,param)
{
	order<-length(param)-1;
	x<-rep(0,length(cp));

	for (k in 1:order) x<-x+param[k]*(cp^(k-1));

	x<-x+param[order+1]*qnorm(cp);
	x
}

normpoly_inv<-function(cp,param)
{
qnormpoly(cp,param)
}


pnormpoly<-function(x,param)
{
# A naive binary search is employed. Error <= 2^(-epsilon).

epsilon<-50;  # max 52

cp<-rep(0.5,length(x));
for (j in seq(-2,-epsilon,by=-1)) cp<-cp+2^j*sign(x-normpoly_inv(cp,param));

cp
}

normpoly_cdf<-function(x,param)
{
pnormpoly(x,param)
}




dnormpoly<-function(x,param)
{
order<-length(param)-1;
cp<-normpoly_cdf(x,param);

denom<-rep(0,length(cp)); 
for (k in 2:order) denom<-denom+param[k]*(k-1)*cp^(k-2);
denom<-denom+param[order+1]*sqrt(2*pi)*exp((qnorm(cp)/sqrt(2))^2); 
density<-1/denom;

# If coefficient for the normal part is zero, the distribution is bounded.
if (param[order+1]==0) density<-density*(x<=sum(param))*(x>=param[1]);

density
}

normpoly_pdf<-function(x,param)
{
dnormpoly(x,param)
}



rnormpoly<-function(n,param)
{
qnormpoly(runif(n),param);
}

normpoly_rnd<-function(n,param)
{
rnormpoly(n,param)
}


covnormpoly4<-function(data)
{
#data is vector
Lcovmat<-Lmomcov(data);
eta2<-1/sqrt(pi);
eta4<-(30*atan(sqrt(2))/pi-9)*eta2;
A<-t(array(c(1,-3,5,3*eta2/eta4,0,6,-30,-6*eta2/eta4,0,0,30,0,0,0,0,1/eta4),c(4,4)))
Acovmat<-A%*%Lcovmat%*%t(A);
Acovmat
}


data2cauchypoly4<-function(data)
{
tlmoms<-t1lmoments(data)
params<-t1lmom2cauchypoly4(tlmoms)
params
}


t1lmom2cauchypoly4<-function(t1lmom)
{
L1<-t1lmom[1];
L2<-t1lmom[2];
L3<-t1lmom[3];
L4<-t1lmom[4];

v1<-0.6978;
v2<-0.3428*v1;

a1<-(25*L4*v1+5*L1*v2-25*L2*v2+63*L3*v2)/(5*v2);
a2<-10*L2-63*L3-10*L4*v1/v2;
a4<-L4/v2; 
a3<-63*L3;
  
param<-c(a1,a2,a3,a4);
param
}


qcauchypoly<-function(cp,param)
{
	order<-length(param)-1;
	x<-rep(0,length(cp));

	for (k in 1:order) x<-x+param[k]*(cp^(k-1));

	x<-x+param[order+1]*tan(pi*(cp-0.5));
	x
}

cauchypoly_inv<-function(cp,param)
{
	qcauchypoly(cp,param)
}


pcauchypoly<-function(x,param)
{
	# A naive binary search is employed. Error <= 2^(-epsilon).

	epsilon<-50;  # max 52

	cp<-rep(0.5,length(x));
	for (j in seq(-2,-epsilon,by=-1)) cp<-cp+2^j*sign(x-cauchypoly_inv(cp,param));

	cp
}

cauchypoly_cdf<-function(x,param)
{
	pcauchypoly(x,param)
}



dcauchypoly<-function(x,param)
{
order<-length(param)-1;
cp<-cauchypoly_cdf(x,param);

denom<-rep(0,length(cp)); 
for (k in 2:order) denom<-denom+param[k]*(k-1)*cp^(k-2);
denom<-denom+param[order+1]*pi*(1/cos(pi*(cp-0.5)))^2;
density<-1/denom;

# If coefficient for the normal part is zero, the distribution is bounded.
if (param[order+1]==0) density<-density*(x<=sum(param))*(x>=param[1]);

density
}

cauchypoly_pdf<-function(x,param)
{
dcauchypoly(x,param)
}


rcauchypoly<-function(n,param)
{
qcauchypoly(runif(n),param);
}

cauchypoly_rnd<-function(n,param)
{
rcauchypoly(n,param)
}

