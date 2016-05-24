ll_gauss.MSAR <-
function(data,theta,gamma,regime) {

T = dim(data)[1]
if(is.null(T)){T=length(data)}
N.samples = dim(data)[2]
if(is.null(N.samples) || is.na(N.samples)){N.samples=1}

d=dim(as.array(data))[3]
if(is.null(d) || is.na(d)){d=1}
data=array(data,c(T,N.samples,d))

order = attributes(theta)$order

gamma = matrix(gamma[,,regime],N.samples,dim(gamma)[2])
A0 = theta$A0[regime,]
A = list()
A[[1]] = theta$A[[regime]][[1]] # order=1
sigma = theta$sigma[[regime]]

f = 0
for(ex in 1:N.samples){
	ar = array(0,c(d,T-order))
	A0.emis = matrix(0,d,T-order)
	for(o in 1:order){
		ar=ar+A[[o]]%*%t(data[((order+1):T)-o,ex,])
	}
	for(i in 1:d){
		A0.emis[i,]=ar[i,]+A0[i]	
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
