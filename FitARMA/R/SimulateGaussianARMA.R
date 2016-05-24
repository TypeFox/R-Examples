"SimulateGaussianARMA" <-
function(phi, theta, n, InnovationVariance=1, UseC=TRUE)
{
    p<-length(phi)
    q<-length(theta)
    r<-max(p,q)
    a<-rnorm(n+q, mean=0, sd=sqrt(InnovationVariance))
    if(p==0 && q==0) return(a)
    z<-numeric(n)
    g<-numeric(0)
    if (p>0) 
        g<-TacvfARMA(phi,theta,p-1)
    if (is.null(g)){ #is null only if non-causal
         warning("Simulating non-stationary stochastic difference equation")
         z[1:p]<-a[1:p]
         } 
    if (p>0)
        if (q>0) {
        psi<-ImpulseCoefficientsARMA(phi,theta,q)[1:q]
        v12<-matrix(rep(0,p*q), nrow=p, ncol=q)
        for (i in 1:p)
            for (j in 1:i) 
                v12[i,j]<-psi[1+j-i]
        id<-matrix(rep(c(1,rep(0,q)),q)[1:q^2], nrow=q, ncol=q)
        v<-cbind(rbind(toeplitz(g), v12),rbind(id,t(v12)))
        sdecomp<-svd(v)
        sqrtv<-(sdecomp$v)%*%diag(sqrt(sdecomp$d))%*%(sdecomp$u)
        za<-crossprod(a[1:(p+q)],sqrtv)
        z[1:p]<-za[1:p]
        z[1:q]<-za[p+(1:q)]
        }
        else 
           z[1:p]<-crossprod(a[q+(1:p)],chol(toeplitz(g)))
    if (UseC){
        if (!is.loaded("GetSimARMA")){
            message("GetSimARMA not loaded. Trying dyn.load ...")
	    return()
            #dyn.load("d:/r/2005/faster/SimGA.dll")
        }
        beta<-c(phi,theta)
        par<-c(n,p,q)
        ans<-.C(GetSimARMA,z,a,beta,par)
        z<-ans[[1]]
    }
    else {
        i<-rep(1:(q+1),n)
        j<-rep(0:(n-1),rep((q+1),n))
        ij<-i+j
        u<-c(crossprod(rev(c(1,-theta)),matrix(a[ij],nrow=(q+1))))
        if(p==0) return(u)
        for (i in (r+1):n) z[i]=u[i]+sum(rev(phi)*z[(i-p):(i-1)])
        }
    z   
}

