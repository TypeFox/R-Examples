ll_gausshn_p0.MSAR <-
function(a,data,gamma,nh.emissions,covar) {
order = 0
T <- dim(as.array(data))[1]
if(is.null(T) || is.na(T)){T <- length(data)}
N.samples <- dim(as.array(data))[2]
if(is.null(N.samples) || is.na(N.samples)){N.samples <- 1}
d <- dim(as.array(data))[3]
if(is.null(d) || is.na(d)){d <- 1}
data <- array(data,c(T,N.samples,d))
ncov=dim(as.array(covar))[3]
if(is.null(ncov) || is.na(ncov)){ncov=1}
f = 0

moy = matrix(a[1:d],d,1)
sigma = exp(a[(d+1):(d+d^2)])
sigma=matrix(sigma,d,d)
#if (length(a)>((order+1)*d^2+d)) {
	par.emis = a[(d^2+d+1):length(a)]
	par.emis = matrix(par.emis,d,ncov)
#	}
#else {par.emis = matrix(0,d,2)}

for (ex in 1:N.samples){
   moyenne = matrix(0,d,T);
   if (d==1) {moyenne[1,]=moyenne[1,]+moy[1]+nh.emissions(as.matrix(covar[,ex,]),par.emis)
   } else { for(i in 1:d){
   	moyenne[i,]=moyenne[i,]+moy[i]+nh.emissions(as.matrix(covar[,ex,]),par.emis)[[i]]
   }}
   # Gaussian density                    
   dn = matrix(pdf.norm(t(data[(1:T),ex,]), moyenne, sigma),T,dim(data)[3])
   dn[dn<1e-15] = 1e-15 ;
   B = log(dn) 
   f = f+sum(gamma[ex,]*t(B));
   }

f=-f;
return(f) 
}
