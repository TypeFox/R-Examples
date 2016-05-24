`SimulateGaussianAR` <-
function(phi, n=100, InnovationVariance=1)
{
    p<-length(phi)
    a<-rnorm(n, mean=0, sd=sqrt(InnovationVariance))
    if(p==0) return(a)
    if (p==1 && phi==1) return(cumsum(a)) #convenient for unit root test
    z<-numeric(n)
    g<-TacvfAR(phi,p-1)
    if (is.null(g)){ #is null only if non-causal
         warning("Simulating non-stationary stochastic difference equation")
         z[1:p]<-a[1:p]
         } 
    if (p>0)
        z[1:p]<-crossprod(a[1:p],chol(toeplitz(g)))
       for (i in (p+1):n) 
            z[i]=a[i]+sum(rev(phi)*z[(i-p):(i-1)])
    z   
}

