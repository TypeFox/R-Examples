diffpci=function(a,b,c,d,n,alpha){
    if (a+b+c+d!=n){return(paste("Caution: a+b+c+d is not equal to n."))}
    if (a+b+c+d==n){

        wald = diffpropci.Wald.mp(c, b, n, 1 - alpha)
        wald.cc = wald.cc(b, c, n, alpha)
        agresti = diffpropci.mp(c, b, n, 1 - alpha)
        tango2 = scoreci.mp(c, b, n, 1 - alpha)
        exact.cond = exact.cond(b, c, n, alpha)
        exact.midp = exact.midp(b, c, n, alpha)
        uncond = uncond(a, b, c, d, n, alpha)
        wilson = wilson(a, b, c, d, n, alpha)
        wilson.cc = wilson.cc(a, b, c, d, n, alpha)
        wilson.phi = wilson.phi(a, b, c, d, n, alpha)
        np.nv = np.nv(a, b, c, d, n, alpha)
        np.t = np.t(a,b,c,d,n,alpha)
                                        #method=c("lower limit","upper limit")
        out=structure(list(Method=c("Wald","Wald.cc","Agresti","Tango", "Exact.cond", "Exact.midp", "Uncond", "Wilson", "Wilson.cc", "Wilson.phi", "np.nv","np.t"),
        estimator=c(wald$estimate,wald.cc$estimate, agresti$estimate,exact.cond$estimate,exact.cond$estimate,exact.midp$estimate,uncond$estimate,wilson$estimate,wilson.cc$estimate,wilson.phi$estimate,np.nv$estimate,np.t$estimate),
        lower=c(wald$conf.int[1],wald.cc$conf.int[1],agresti$conf.int[1], tango2$conf.int[1],exact.cond$conf.int[1],exact.midp$conf.int[1],uncond$conf.int[1],wilson$conf.int[1],wilson.cc$conf.int[1],wilson.phi$conf.int[1],np.nv$conf.int[1],np.t$conf.int[1]),
        upper=c(wald$conf.int[2],wald.cc$conf.int[2],agresti$conf.int[2], tango2$conf.int[2],exact.cond$conf.int[2],exact.midp$conf.int[2],uncond$conf.int[2],wilson$conf.int[2],wilson.cc$conf.int[2],wilson.phi$conf.int[2],np.nv$conf.int[2],np.t$conf.int[2])),
        .Names=c("Method","Estimator","lower limit", "upper limit"),
        row.names = c(NA,12L), class = "data.frame")
        return(out)
    }
}
