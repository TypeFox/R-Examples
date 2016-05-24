elltheta <-
function(theta, p0, r, tol=1e-7, maxit=500){
    n<-length(r)
    del<-1
    it<-0
    p<-p0
    while(del>tol && it<maxit){
        res<-H.Binv(theta,p, r)
        p.nu<-p-res$Binv%*%res$H
        del<-sum(abs(res$H))+sum(abs(p-p.nu))
        p<-p.nu
        #cat(it, del,  "theta=", theta, "sum of p=", sum(p), "\n")
        it<-it+1
    }
    ell<-res$ell
    Fisher<-res$Fisher
    list(ell=ell,  p=p)
}    
