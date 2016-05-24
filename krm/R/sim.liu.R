# Liu et al 2008
sim.liu.2008=function (n, a, seed=NULL) {
    
    if(!is.null(seed)) set.seed(seed)
    
    requireNamespace("MASS")
    z=MASS::mvrnorm(n=n, mu=rep(0,5), Sigma=diag(rep(1,5)))
    h.z = 2*(z[,1]-z[,2])^2 + z[,2]*z[,3] + 3*sin(2*z[,3])*z[,4] + z[,5]^2 + 2*cos(z[,4])*z[,5]    
    x=z[,1]+rnorm(n, mean=0)/2 # version 3 and before, mean is 0.5
    
    y=rbern(n,expit(x + a*h.z))
    
    data.frame(y,x,z=z) # "y"   "x"   "z.1" "z.2" "z.3" "z.4" "z.5"
    
}



# Liu et al 2007
sim.liu.2007=function (n, a, seed=NULL) {
    
    if(!is.null(seed)) set.seed(seed)
    
    requireNamespace("MASS")
    z=MASS::mvrnorm(n=n, mu=rep(0,5), Sigma=diag(rep(1,5)))
    h.z = 2*cos(z[,1]) - 3*z[,2]^2 + 2*exp(-z[,3])*z[,4] - 1.6*sin(z[,5])*cos(z[,3]) + 4*z[,1]*z[,5]
    x=3*cos(z[,1])+2*rnorm(n, mean=0) 
    
    y=x + a*h.z + rnorm(n)
    
    data.frame(y,x,z=z) # "y"   "x"   "z.1" "z.2" "z.3" "z.4" "z.5"
    
}
#sim.liu.2007 (n=5, a=1, seed=1) 
