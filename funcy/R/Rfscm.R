#
#  Copyright (C) 2013  Huijing Jiang & Nicoleta Serban
#  
#

Rfscm=function(y,nb,phi,psi,z,max_iter){

    ##Initialize beta##
    n=nrow(y)
    m=nrow(y)
    p=ncol(phi)
    q=ncol(psi)
    gamma=rep(0,q)
    beta=estimata.beta(y,z,phi,psi,gamma)

    ##Initialize gamma##
    sigma_e=0
    sigma_s=1
    gamma=impute.gamma(y,z,phi,psi,beta,sigma_e,sigma_s)

    ##Initialize variance components
    sigma=estimate.var(y,z,phi,psi,beta,gamma,sigma_e,sigma_s)
    sigma_e=sigma$sigma_e
    #cat("Sigma_epsilon is", sigma_e, "\n")
    sigma_s=sigma$sigma_s
    #cat("Sigma_s is", sigma_s, "\n")
    Vgamma=sigma$Vgamma    
    
    ##Initialize pi
    theta=optimize(log.prob_z_y,c(0, 1),z,nb)$minimum
    #cat("Gibbs parameter is ", theta, "\n")
    prob_z=gibbs(theta,z,nb)

    iter=1
    max_iter_mc=100

    while(iter<=max_iter){

        #cat("The", iter,"th iteration:\n")
        ##E-step#
        gamma=impute.gamma(y,z,phi,psi,beta,sigma_e,sigma_s)
        z=impute.z(y,z,phi,psi,beta,gamma,prob_z,sigma_e,Vgamma,100)
        k <- max(z)
        if(max(z)!= length(unique(z))){
            warning("class number was reduced!")
            res <- label2lowerk(z)
            z <- res$cluster
            k <- res$k
        }
           
        ##M-step#
        beta=estimata.beta(y,z,phi,psi,gamma)

        sigma=estimate.var(y,z,phi,psi,beta,gamma,sigma_e,sigma_s)
        sigma_e=sigma$sigma_e
        #cat("Sigma_epsilon is", sigma_e, "\n")
        sigma_s=sigma$sigma_s
        #cat("Sigma_s is", sigma_s, "\n")
        Vgamma=sigma$Vgamma    
        
        theta=optimize(log.prob_z_y,c(0, 1),z,nb)$minimum
        #cat("Gibbs parameter is ", theta, "\n")
        prob_z=gibbs(theta,z,nb)

        ##update iteration
        iter=iter+1
    }

    return(resulst=list(z=z, k=k, beta=beta,gamma=gamma,sigma_e=sigma_e,sigma_s=sigma_s,theta=theta, prob_z=prob_z))

}




impute.z=function(y,z,phi,psi,beta,gamma,pi,sigma_e,Vgamma,max_iter_mc){
    
    n=nrow(y)
    m=ncol(y)
    c=ncol(beta)
    
    ##yres=y
    pi_y=matrix(0,n,c)
    
    r=0
    max_iter_mc=100 
    while(r<max_iter_mc){

        pi_y_r=matrix(0,n,c)
        gamma_r=mvrnorm(1,gamma,Vgamma)
        tau=psi%*%gamma_r
        
        for(s_j in 1:n){
            ##yres[s_j,]=yres[s_j,]-phi%*%alpha
            for(c_i in 1:c){
                prob_y_c=exp(-0.5/sigma_e*sum((y[s_j,]-tau[s_j]-phi%*%beta[,c_i])^2))
                if(prob_y_c==0)
                    prob_y_c <- 0.00001
                pi_y_r[s_j,c_i]=pi[s_j,c_i]*prob_y_c
            }##end of c            
            
        }##end of s
        
        pi_y_r=t(apply(pi_y_r,1,function(x) x/sum(x)))
        pi_y=pi_y+pi_y_r
        r=r+1
        
    }##

    pi_y=pi_y/max_iter_mc
    ##apply hard clustering
    
    z=apply(pi_y,1,function(x) which(x==max(x)))
    
    return(z)
    
}

log.prob_z_y=function(theta,z,nb){

    prob_z_y=0 
    n=length(z)
    c <- max(z)
    
    for(s_j in 1:n){
        log_pi=rep(0,c)
        n_c=table(z[nb[s_j,]])               
        log_pi[as.numeric(names(n_c))]=theta*n_c-log(sum(exp(theta*n_c)))  
        prob_z_y=prob_z_y-log_pi[z[s_j]]
    }
    return(prob_z_y)

}

gibbs=function(theta,z,nb){

    n=length(z)
    uniz <- sort(unique(z))
    c=length(uniz)
    pi=matrix(0,n,c)   
    colnames(pi) <- uniz
    
    for(s_j in 1:n){
        n_c=table(z[nb[s_j,]])
        pi[s_j,names(n_c)]=exp(theta*n_c)
    }
    
    pi=apply(pi,1,function(x) x/sum(x))

    return(t(pi))
    
}

estimata.beta=function(y,z,phi,psi,gamma){

    c <- max(z)
    n=nrow(y)
    m=ncol(y)
    p=ncol(phi)
    q=ncol(psi)
    n_k=table(z)
    lambda=rep(0,p)
    
    ##remove spatial effects
    tau=psi%*%gamma
    yres=y
    for(j in 1:n)
        yres[j,]=yres[j,]-tau[j]
    
    for(j in 1:n)
        lambda=lambda+t(phi)%*%yres[j,]/n_k[names(n_k)==z[j]]

    lambda=lambda/sum(1/n_k)
    beta=matrix(rep(0,c*p),p,c)

    for(j in 1:n)
        beta[,z[j]]=beta[,z[j]]+t(phi)%*%yres[j,]
    
    beta=apply(beta,2,function(x) x-lambda)
    for(c_i in 1:c)
        beta[,c_i]=solve(t(phi)%*%phi)%*%beta[,c_i]/n_k[c_i]
    
    return(beta)
    
}

impute.gamma=function(y,z,phi,psi,beta,sigma_e,sigma_s){

    n=nrow(y)
    m=ncol(y)
    p=ncol(phi)
    q=ncol(psi)

    V=m*t(psi)%*%psi
    diag(V)=diag(V)+sigma_e/sigma_s
    Vinv=solve(V)
    
    yres=y
    ##remove cluster effects
    for(j in 1:n)
        yres[j,]=yres[j,]-phi%*%beta[,z[j]]
    
                                        #Vgamma=sigma_e*Vinv
    gamma=Vinv%*%t(psi)%*%yres

    gamma=apply(gamma,1,sum)
    return(gamma)

}

estimate.var=function(y, z, phi, psi, beta, gamma, sigma_e, sigma_s){
    
    n=nrow(y)
    m=ncol(y)
    p=ncol(phi)
    q=ncol(psi)
    
    ##estimate gamma_s
    V=m*t(psi)%*%psi
    diag(V)=diag(V)+sigma_e/sigma_s
    Vinv=solve(V)
    
    Vgamma=sigma_e*Vinv
    
    sigma_s=mean(gamma^2)
    sigma_s=sigma_s+mean(diag(Vgamma))

    yres=y
    tau=psi%*%gamma

    for(j in 1:n){
        yres[j,]=yres[j,]-phi%*%beta[,z[j]]
        yres[j,]=yres[j,]-tau[j]
    }

    sigma_e=mean(yres^2)
    #cat(sigma_e,"\n")
    sigma_e=sigma_e+mean(diag(psi%*%Vgamma%*%t(psi)))
    #cat(sigma_e,"\n")
    
    return(sigma=list(sigma_e=sigma_e,sigma_s=sigma_s,Vgamma=Vgamma))

}
