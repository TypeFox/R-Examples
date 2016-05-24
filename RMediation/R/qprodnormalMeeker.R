qprodnormalMeeker <-
function(p, mu.x, mu.y, se.x, se.y, rho=0, lower.tail=TRUE){
    max.iter=1000
    eps <- 1; #intial value for precision
    mu.a <- mu.x/se.x #Rescaling distribution of a so that SD of a is 1
    mu.b <- mu.y/se.y #Rescaling distribution of b so that SD of b is 1
    se.ab<-sqrt(1+mu.a^2+mu.b^2+2*mu.a*mu.b*rho+rho^2)
    s.a.on.b <- sqrt(1-rho^2) #SD of b conditional on a
    if (lower.tail==FALSE) {
        u0<- mu.a*mu.b +6*se.ab #upper
        l0=mu.a*mu.b- 6*se.ab
        alpha<- 1- p
    }
    else{
        l0<-mu.a*mu.b - 6*se.ab #lower
        u0<-mu.a*mu.b + 6*se.ab
        alpha<- p
    }
    gx=function(x, z) {
        mu.a.on.b <- mu.a+rho*(x-mu.b) #mean of a conditional on b
        integ <- pnorm(sign(x)*(z / x-mu.a.on.b)/s.a.on.b)*dnorm(x-mu.b)
        return(integ)
    }
    fx<-function(z){
    	return(integrate(gx,lower=-Inf,upper=Inf, z=z)$value-alpha)
    }
    p.l<-fx(l0)
    p.u<-fx(u0)
    iter=0
    while (p.l>0) {
        iter=iter+1
        l0=l0-0.5*se.ab
        p.l<-fx(l0)
        if (iter> max.iter )
            {
                cat(" No initial valid lower bound interval!\n")
                return(list(q=NA,error=NA))
            }
    }
    iter=0 #Reset iteration counter
    while (p.u<0) {
        iter=iter+1
        u0=u0+0.5*se.ab
        p.u<-fx(u0)
        if (iter> max.iter ){
            cat(" No initial valid upper bound interval!\n")
            return(list(q=NA,error=NA))
        }
   }
    res<-uniroot(fx,c(l0,u0))
    new<-res$root
    new <- new*se.x*se.y
    error.new<-res$estim.prec
    return(list(q=new,error=error.new))   #returns quantile corresponding to p
}
