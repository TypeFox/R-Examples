solve.beta <-
function(beta.ini,Rt.ini,obs.t,delta,z,wt,r,dx,iter.max){
    zc = t(t(z) - apply(z,2,mean))
    ix = order(obs.t)
    ord.delta = delta[ix]
    ord.z = as.matrix(zc[ix,])
    colnames(ord.z)= colnames(z)
    ord.wt = wt[ix]
    p = length(beta.ini)
    n = length(Rt.ini)

    dif = 1
    iter = 0
    while(dif > dx & iter <= iter.max){
        ord.bz = as.numeric(ord.z%*%beta.ini)  
        s0 = S0.fun(ord.bz,Rt.ini,ord.wt,r,n)
        Rt = Rt.fun(ord.delta,ord.wt,s0)

        Ubeta = Ubeta.fun(ord.delta,ord.z,ord.bz,Rt,ord.wt,r)
        Hbeta = Hbeta.fun(ord.z,ord.bz,Rt,ord.wt,r,p)        

        beta = beta.ini + as.numeric(solve(Hbeta)%*%Ubeta)

        dif = max(abs(beta-beta.ini))
        iter = iter + 1
        beta.ini = beta
        Rt.ini = Rt
    }
    converged = as.numeric(dif > dx)

    return(list(Beta=beta,converged=converged,iter,data=data.frame(ord.time=sort(obs.t),Rt=Rt)))
}