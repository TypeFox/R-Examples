confid.int.theta <-
function(x, y, method=c("chi-sq", "simulate"), conf.level = 0.95, grd = .001,  
                            B = 1000, tol=1e-7, maxit=500){
    method <- match.arg(method)
    alfa<-1-conf.level
    m<-length(x)
    n<-length(y)
    N<-m+n
    lambda<-m/N
    z<-c(y,x)
    r<-rank(z, ties.method = "max")[1:n]
    theta0<-dd.est(x,y)
    p0<-phi(N, theta0, lambda)/N
    res<-mrle.sporm(x, y, theta0, p0)
    phat<-res$p
    theta.hat<-res$theta
    ell<-res$ell
    if(method == "chi-sq"){
        Calfa<-qchisq(conf.level,1)
    }
    else{
        LR.b<-NULL
        sim.c<-function(){
            u<-runif(m)
            v<-runif(n)
            v<-v/(theta.hat-(theta.hat-1)*v)
            rv<-rank(c(v,u), ties.method = "max")[1:n]
            theta0<-dd.est(u,v)
            p0<-phi(N, theta0, lambda)/N
            ell.star<-try(mrle.sporm(u, v, theta0, p0)$ell, TRUE) 
            p0<-phi(N, theta.hat, lambda)/N
            2*(ell.star-try(elltheta(theta.hat, p0, rv, tol, maxit)$ell, TRUE))
        }
        LR.b <- lapply(1:B, function(i) try(sim.c(), TRUE))
        LR.b <-unlist(LR.b[sapply(LR.b, function(x) !inherits(x, "try-error"))])
        Calfa<-quantile(LR.b, conf.level)
    }
    theta<-theta.hat+grd
    p<-phat
    LR.theta<-2*(ell-elltheta(theta, p, r, tol, maxit)$ell)
    while(LR.theta<Calfa){
        theta<-theta+grd
        tmp<-elltheta(theta, p, r, tol, maxit)
        LR.theta<-2*(ell-tmp$ell)
        p<-tmp$p
    }
    theta.U<-theta 
    theta<-theta.hat-grd
    LR.theta<-2*(ell-elltheta(theta, p, r, tol, maxit)$ell)
    while(LR.theta<Calfa & theta>0){
        theta<-theta-grd
        tmp<-elltheta(theta, p, r, tol, maxit)
        LR.theta<-2*(ell-tmp$ell)
        p<-tmp$p
    }
    theta.L<-theta 
    list(theta.L=theta.L, theta.U=theta.U, theta.hat=theta.hat, Calpha=Calfa)
}
