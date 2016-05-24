ll_gausshn.MSAR <-
function(a,data,gamma,order,nh.emissions,covar) {

T = dim(data)[1]
if(is.null(T)){T=length(data)}
N.samples = dim(data)[2]
if(is.null(N.samples) || is.na(N.samples)){N.samples=1}

d=dim(as.array(data))[3]
if(is.null(d) || is.na(d)){d=1}
ncov=dim(as.array(covar))[3]
if(is.null(ncov) || is.na(ncov)){ncov=1}
data=array(data,c(T,N.samples,d))

f = 0

Atmp = as.vector(a[(1:(order*d^2))])
A=NULL
for(o in 1:order){
	A[[o]]=matrix(Atmp[((o-1)*d^2+1):(o*d^2)],d,d)
}

A0 = a[(order*d^2+1):(order*d^2+d)]

sigmatmp =exp(a[(order*d^2+d+1):((order+1)*d^2+d)])
#sigmatmp =(a[(order*d^2+d+1):((order+1)*d^2+d)])
#write(sigmatmp, "", sep = "\t")
#write(A0, "", sep = "\t")

sigma=matrix(sigmatmp,d,d)

#if (length(a)>((order+1)*d^2+d)) {
	par.emis = a[((order+1)*d^2+d+1):length(a)]
#	par.emis = matrix(par.emis,d,2)
	par.emis = matrix(par.emis,d,ncov)
#	}
#par.emis[2] = par.emis[2] %% (2*pi)
#else {par_emis = matrix(0,d,2)}

for(ex in 1:N.samples){
	ar = array(0,c(d,T-order))
	A0.emis = matrix(0,d,T-order)
	for(o in 1:order){
		ar=ar+A[[o]]%*%t(data[((order+1):T)-o,ex,])
	}
	for(i in 1:d){
		femis = nh.emissions(as.matrix(covar[(order + 1):T, ex,]),par.emis)
		A0.emis[i,]=ar[i,]+A0[i]+femis[i,]	
	}
# Gaussian pdf
dn=pdf.norm(t(data[(order+1):T,ex,]),A0.emis,sigma)
dn[dn<1e-10] = 1e-10
B = log(dn) 
f = f+sum(gamma[ex,]*B[1:length(gamma[ex,])])
}

f=-f;
return(f) 
}
