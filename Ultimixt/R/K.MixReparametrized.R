K.MixReparametrized <-
function(xobs, k, alpha0, alpha, Nsim){
    alpha0=rep(alpha0, k)
    alpha1=alpha
    n=length(xobs)
    p0=rep(1/k, k)
    ta2=rep(1/(k+1),k)
    tau0=sqrt(ta2)
    fhi0=sqrt(1/(k+1))
    thita0=1
    nu0=1
    
    sd_mu=.1
    sd_sig=.05
    eps1=500
    eps2=500
    eps=.1
    eps_nu=.1
    betta=rep(1, 4)
    T=30000
    m=length(betta)
    naccept1=rep(0,m)
    naccept2=rep(0,m)
    naccept3=rep(0,m)
    naccept4=rep(0,m)
    naccept5=rep(0,m)
    naccept6=rep(0,m)
    na1=0
    na2=0
    na3=0
    na4=0
    na5=0
    na6=0
    adapt_mu_ra=rep(0,T)
    adapt_sd_mu=rep(0,T)
    adapt_sig_ra=rep(0,T)
    adapt_sd_sig=rep(0,T)
    adapt_p_ra=rep(0,T)
    adapt_eps_p=rep(0,T)
    adapt_degata_ra=rep(0,T)
    adapt_eps_degata=rep(0,T)
    adapt_thita_ra=rep(0,T)
    adapt_eps_thita=rep(0,T)
    adapt_nu_ra=rep(0,T)
    adapt_eps_nu=rep(0,T)
    
    n=length(xobs)
    meanobs=mean(xobs)
    
    mui=rep(mean(xobs), m)
    output_mu_s=matrix(mean(xobs), nrow=T, ncol=m)
    
    ss2=rep(var(xobs), m)
    output_sigma2_s=matrix(var(xobs), nrow=T, ncol=m)
    ss=rep(sd(xobs), m)
    output_sigma_s=sqrt(output_sigma2_s)
    
    pp=matrix(p0, nrow=m, ncol=k)
    output_p_s=array(p0, dim=c(T,m,k))
    
    post_dist_mu=matrix(0, nrow=T, ncol=m)
    post_dist_sigma=matrix(0, nrow=T, ncol=m)
    post_dist_p=matrix(0, nrow=T, ncol=m)
    post_dist_fhi=matrix(0, nrow=T, ncol=m)
    post_dist_thita=matrix(0, nrow=T, ncol=m)
    post_dist_nu=matrix(0, nrow=T, ncol=m)
    
    # if k>2
    if(k>2){
        thi=matrix(thita0, nrow=m, ncol=(k-2))
        output_thita_s=array(thita0, dim=c(T,m,(k-2)))
    }else{output_thita_s=0}
    
    nu=matrix(nu0, nrow=m, ncol=(k-1))
    output_nu_s=array(nu0, dim=c(T,m,(k-1)))
    
    fhi=rep(fhi0, m)
    output_fhi_s=matrix(fhi0, nrow=T, ncol=m)
    
    component_mean=array(0, dim=c(T,m,k))
    component_sigma=array(0, dim=c(T,m,k))
    F=array(0, dim=c(m, k, (k-1)))
    
    a_adapt=0 #batch number of 100 iterations
    a=1
    
    den_dirichlet<-function (x, alpha)
    {
        dirichlet1 <- function(x, alpha) {
            logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
            s <- sum((alpha - 1) * log(x))
            exp(sum(s) - logD)
        }
        if (!is.matrix(x))
        if (is.data.frame(x))
        x <- as.matrix(x)
        else x <- t(x)
        if (!is.matrix(alpha))
        alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x),
        byrow = TRUE)
        if (any(dim(x) != dim(alpha)))
        stop("Mismatch between dimensions of x and alpha in ddirichlet().\n")
        pd <- vector(length = nrow(x))
        for (i in 1:nrow(x)) pd[i] <- dirichlet1(x[i, ], alpha[i,
        ])
        pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
        pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0
        return(pd)
    }
    rdirichlet=function(n, alpha){
        l <- length(alpha)
        x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
        sm <- x %*% rep(1, l)
        return(x/as.vector(sm))
    }
    ############## Orthonormal basis:
    vector_gamma=function(p, theta){
        F[,1,1]=-sqrt(p[,2])/sqrt(p[,1]+p[,2])
        F[,2,1]=sqrt(p[,1])/sqrt(p[,1]+p[,2])
        if(k>2){
            for(s in 2:(k-1)){
                sp=sqrt(rowSums(p[,1:s]))
                st=p[,1:s]*p[,(s+1)]
                s1=-sqrt(st)/sp
                
                F[,1:s,s]=s1
                F[,(s+1),s]=sp
                
                sF2=matrix(F[,,s]^2, ncol=k)
                sF=rowSums(sF2)
                ssF=sqrt(sF)
                
                F[,,s]=sweep(matrix(F[,,s], ncol=k),MARGIN=1,ssF,`/`)
            }
            ############## angles:
            sin_theta=matrix(sin(theta), ncol=(k-2))
            cos_theta=matrix(cos(theta), ncol=(k-2))
            
            M_theta=array(1, dim=c((k-1), (k-2),m))
            v=1
            repeat{
                M_theta[,1,v]=c(cos_theta[v,], sin_theta[v,(k-2)])
                # if k>3:
                if(k>3){
                    for(l in 2:(k-2)){
                        M_theta[l,2:l,v]=sin_theta[v,1:(l-1)]
                    }
                    M_theta[(k-1),2:(k-2),v]=sin_theta[v,1:(k-3)]
                }
                
                if(v==m)
                break
                v=v+1
            }
            
            vec_angles=matrix(0, ncol=(k-1), nrow=m)
            v=1
            repeat{
                if(k>3){
                    vec_angles[v,]=apply(M_theta[,,v], 1, prod)
                }else{
                    vec_angles[v,]=M_theta[,,v]
                }
                
                if(v==m)
                break
                v=v+1
            }
            ############## gammas:
            vec_gamma=matrix(0, ncol=k, nrow=m)
            M_gamma=array(0, dim=c(k, (k-1), m))
            v=1
            repeat{
                M_gamma[,,v]=sweep(F[v,,],MARGIN=2,vec_angles[v,],`*`)
                vec_gamma[v,]=rowSums(M_gamma[,,v])
                if(v==m)
                break
                v=v+1
            }
            
        }else{
            vec_gamma=matrix(0, ncol=k, nrow=m)
            vec_gamma[,1]=F[,1,1]
            vec_gamma[,2]=F[,2,1]
        }
        return(vec_gamma)
    }
    
    ############## angles 2:
    vector_angles_nu=function(nu){
        
        sin_nu=matrix(sin(nu), ncol=(k-1))
        cos_nu=matrix(cos(nu), ncol=(k-1))
        
        M_nu=array(1, dim=c(k, (k-1),m))
        v=1
        repeat{
            M_nu[,1,v]=c(cos_nu[v,], sin_nu[v,(k-1)])
            # if k>2:
            if(k>2){
                for(l in 2:(k-1)){
                    M_nu[l,2:l,v]=sin_nu[v,1:(l-1)]
                }
                M_nu[k,2:(k-1),v]=sin_nu[v,1:(k-2)]
            }
            
            if(v==m)
            break
            v=v+1
        }
        
        vec_angles_nu=matrix(0, ncol=k, nrow=m)
        v=1
        repeat{
            if(k>2){
                vec_angles_nu[v,]=apply(M_nu[,,v], 1, prod)
            }else{
                vec_angles_nu[v,]=M_nu[,,v]
            }
            
            if(v==m)
            break
            v=v+1
        }
        return(vec_angles_nu)
    }
    
    rept=1
    adapt_condition_mu=TRUE
    adapt_condition_sigma=TRUE
    adapt_condition_p=TRUE
    adapt_condition_fhi=TRUE
    adapt_condition_thita=TRUE
    adapt_condition_nu=TRUE
    adapt_condition=TRUE
    T1=1
    T_stop=1000
    # T_stop=15000
    criter_condition=FALSE
    a_mu=0
    a_sigma=0
    a_p=0
    a_fh=0
    a_theta=0
    a_nu=0
    repeat{
        
        #################################
        for(t in 1:(T-1)){
            
            T1[rept==1]=T1[rept==1]+1
            fhi_gamma=output_fhi_s[t,]*vector_gamma(output_p_s[t,,], output_thita_s[t,,])
            ############## etas:
            fhi_eta=sqrt(1-output_fhi_s[t,]^2)*(vector_angles_nu(output_nu_s[t,,]))
            ################# Hasting independent proposal for mu:
            comp_mean=output_mu_s[t,]+(sweep(fhi_gamma,MARGIN=1,output_sigma_s[t,],`*`)/sqrt(output_p_s[t,,]))
            comp_sigma=sweep(fhi_eta,MARGIN=1,output_sigma_s[t,],`*`)/sqrt(output_p_s[t,,])
            component_mean[t,,]=comp_mean
            component_sigma[t,,]=comp_sigma
            v=1
            mix_prob=matrix(0,ncol=m, nrow=n)
            repeat{
                mixprob=apply(cbind(comp_mean[v,], comp_sigma[v,], output_p_s[t,,][v,]), 1, function(u)u[3]*dnorm(xobs, u[1], u[2]))
                mixprob=rowSums(mixprob)
                mix_prob[,v]=mixprob
                if(v==m)
                break
                v=v+1
            }
            cur=mix_prob
            lcur=log(cur)
            slcur= colSums(lcur, na.rm=FALSE, dims=1)
            curlik=slcur
            if(adapt_condition==TRUE){
                prop1=output_mu_s[t,]+rnorm(m, 0,sd_mu)
                log_mix_mut=dnorm(output_mu_s[t,], prop1, sd_mu, log=TRUE)
                log_mix_mu_prop=dnorm(prop1, output_mu_s[t,], sd_mu, log=TRUE)}else{
                    prop1=mean(xobs)+rnorm(m, 0,sd_mu)
                    meanobs=mean(xobs)
                    log_mix_mut=dnorm(output_mu_s[t,], meanobs, sd_mu, log=TRUE)
                    log_mix_mu_prop=dnorm(prop1, meanobs, sd_mu, log=TRUE)
                }
                comp_mean=prop1+(sweep(fhi_gamma,MARGIN=1,output_sigma_s[t,],`*`)/sqrt(output_p_s[t,,]))
                v=1
                mix_prob=matrix(0,ncol=m, nrow=n)
                repeat{
                    mixprob=apply(cbind(comp_mean[v,], comp_sigma[v,], output_p_s[t,,][v,]), 1, function(u)u[3]*dnorm(xobs, u[1], u[2]))
                    mixprob=rowSums(mixprob)
                    mix_prob[,v]=mixprob
                    if(v==m)
                    break
                    v=v+1
                }
                propl=mix_prob
                lpropl=log(propl)
                slpropl=colSums(lpropl, na.rm=FALSE, dims=1)
                proplik=slpropl
                
                log_ro_mu=betta*proplik-betta*curlik+log_mix_mut-log_mix_mu_prop
                Accept_mu=(log(runif(length(betta)))<(log_ro_mu))
                
                mui[Accept_mu]=prop1[Accept_mu]
                output_mu_s[t+1,]=mui
                
                slcur[Accept_mu]=slpropl[Accept_mu]
                post_dist_mu[t,]=exp(slcur)
                
                naccept1[Accept_mu]=naccept1[Accept_mu]+1
                na1[Accept_mu[1]]=na1[Accept_mu[1]]+1
                a_adapt[a_adapt<50]=a_adapt[a_adapt<50]+1
                
                if(a>300){men_adapt_mu=median(adapt_mu_ra[(a-10):a])
                    cond=(men_adapt_mu>.35&men_adapt_mu<.55)
                    adapt_condition_mu[cond==TRUE]=FALSE
                }
                if(adapt_condition_mu==TRUE){
                    ac_mu_ra=na1/50
                    adapt_mu_ra[a][a_adapt==50]=ac_mu_ra
                    ac_cond_mu1=(ac_mu_ra<.44)
                    ac_cond_mu2=(ac_mu_ra>.44)
                    ls_mu=log(sd_mu[1])
                    ls_mu[ac_cond_mu1&a_adapt==50]=ls_mu[ac_cond_mu1&a_adapt==50]-min(.01,1/sqrt(t))
                    ls_mu[ac_cond_mu2&a_adapt==50]=ls_mu[ac_cond_mu2&a_adapt==50]+min(.01,1/sqrt(t))
                    sd_mu[1]=exp(ls_mu)
                    adapt_sd_mu[a][a_adapt==50]=sd_mu[1]
                    a_mu[a_adapt==50]=a_mu[a_adapt==50]+1}else{adapt_mu_ra[a][a_adapt==50]=adapt_mu_ra[a-1][a_adapt==50]
                        adapt_sd_mu[a][a_adapt==50]=adapt_sd_mu[a-1][a_adapt==50]}
                    ################# Hasting independent proposal for sigma:
                    curlik=slcur-log(output_sigma_s[t,])
                    if(adapt_condition==TRUE){
                        sig=exp(log(output_sigma_s[t,])+rnorm(m,0, sd_sig))
                        log_mix_sigt=-log(output_sigma_s[t,])+dnorm(log(output_sigma_s[t,]),log(sig),sd_sig, log=TRUE)
                        log_mix_sig_prop=-log(sig)+dnorm(log(sig),log(output_sigma_s[t,]),sd_sig, log=TRUE)
                    }else{
                        sig=exp(log(sd(xobs))+rnorm(m,0, sd_sig))
                        sdobs=sd(xobs)
                        log_mix_sigt=-log(output_sigma_s[t,])+dnorm(log(output_sigma_s[t,]),log(sdobs),sd_sig, log=TRUE)
                        log_mix_sig_prop=-log(sig)+dnorm(log(sig),log(sdobs),sd_sig, log=TRUE)
                    }
                    prop2=sig^2
                    comp_mean=output_mu_s[t+1,]+(sweep(fhi_gamma,MARGIN=1,sig,`*`)/sqrt(output_p_s[t,,]))
                    comp_sigma=sweep(fhi_eta,MARGIN=1,sig,`*`)/sqrt(output_p_s[t,,])
                    
                    v=1
                    mix_prob=matrix(0,ncol=m, nrow=n)
                    repeat{
                        mixprob=apply(cbind(comp_mean[v,], comp_sigma[v,], output_p_s[t,,][v,]), 1, function(u)u[3]*dnorm(xobs, u[1], u[2]))
                        mixprob=rowSums(mixprob)
                        mix_prob[,v]=mixprob
                        if(v==m)
                        break
                        v=v+1
                    }
                    
                    propl=mix_prob
                    lpropl=log(propl)
                    slpropl=colSums(lpropl, na.rm=FALSE, dims=1)
                    proplik=slpropl-log(sig)
                    
                    log_ro_sigma=betta*proplik-betta*curlik-log_mix_sig_prop+log_mix_sigt
                    
                    Accept_sig=(log(runif(length(betta)))<(log_ro_sigma))
                    
                    ss2[Accept_sig]=prop2[Accept_sig]
                    ss[Accept_sig]=sig[Accept_sig]
                    output_sigma2_s[t+1,]=ss2
                    output_sigma_s[t+1,]=ss
                    
                    slcur[Accept_sig]=slpropl[Accept_sig]
                    post_dist_sigma[t,]=exp(slcur-log(output_sigma_s[t+1,]))
                    naccept2[Accept_sig]=naccept2[Accept_sig]+1
                    na2[Accept_sig[1]]=na2[Accept_sig[1]]+1
                    if(a>300){men_adapt_sigma=median(adapt_sig_ra[(a-10):a])
                        cond=(men_adapt_sigma>.35&men_adapt_sigma<.55)
                        adapt_condition_sigma[cond==TRUE]=FALSE
                    }
                    if(adapt_condition_sigma==TRUE){
                        ac_sig_ra=na2/50
                        adapt_sig_ra[a][a_adapt==50]=ac_sig_ra
                        ac_cond_sig1=(ac_sig_ra<.44)
                        ac_cond_sig2=(ac_sig_ra>.44)
                        ls_sig=log(sd_sig[1])
                        ls_sig[ac_cond_sig1&a_adapt==50]=ls_sig[ac_cond_sig1&a_adapt==50]-min(.01,1/sqrt(t))
                        ls_sig[ac_cond_sig2&a_adapt==50]=ls_sig[ac_cond_sig2&a_adapt==50]+min(.01,1/sqrt(t))
                        sd_sig[1]=exp(ls_sig)
                        adapt_sd_sig[a][a_adapt==50]=sd_sig[1]
                        a_sigma[a_adapt==50]=a_sigma[a_adapt==50]+1}else{adapt_sig_ra[a][a_adapt==50]=adapt_sig_ra[a-1][a_adapt==50]
                            adapt_sd_sig[a][a_adapt==50]=adapt_sd_sig[a-1][a_adapt==50]}
                        ######################### Hasting independent proposal for nus:
                        curlik=slcur
                        
                        pido=pi/2
                        nu_prop=matrix(0, nrow=m, ncol=(k-1))
                        if(k>2){number=m*(k-2)
                            nu_prop[,1:(k-2)]=matrix(runif(number,0,pido), nrow=m, ncol=(k-2))}
                        nu_prop[,(k-1)]=runif(m, 0, pido)
                        ############## etas:
                        fhi_eta=sqrt(1-output_fhi_s[t,]^2)*(vector_angles_nu(nu_prop))
                        ################# Hasting independent:
                        comp_mean=output_mu_s[t+1,]+(sweep(fhi_gamma,MARGIN=1,output_sigma_s[t+1,],`*`)/sqrt(output_p_s[t,,]))
                        comp_sigma=sweep(fhi_eta,MARGIN=1,output_sigma_s[t+1,],`*`)/sqrt(output_p_s[t,,])
                        
                        v=1
                        mix_prob=matrix(0,ncol=m, nrow=n)
                        repeat{
                            mixprob=apply(cbind(comp_mean[v,], comp_sigma[v,], output_p_s[t,,][v,]), 1, function(u)u[3]*dnorm(xobs, u[1], u[2]))
                            mixprob=rowSums(mixprob)
                            mix_prob[,v]=mixprob
                            if(v==m)
                            break
                            v=v+1
                        }
                        
                        propl=mix_prob
                        lpropl=log(propl)
                        slpropl=colSums(lpropl, na.rm=FALSE, dims=1)
                        proplik=slpropl
                        
                        log_ro_nu=betta*proplik-betta*curlik
                        Accept_nu=(log(runif(length(betta)))<(log_ro_nu))
                        
                        nu[Accept_nu,]=nu_prop[Accept_nu,]
                        output_nu_s[t+1,,]=nu
                        
                        slcur[Accept_nu]=slpropl[Accept_nu]
                        
                        if(k>2){
                            ######################### Hasting independent proposal for theta:
                            curlik=slcur
                            
                            dopi=2*pi
                            thita=matrix(0, nrow=m, ncol=(k-2))
                            if(k>3){
                                number=m*(k-3)
                                thita[,1:(k-3)]=matrix(runif(number,0,pi), nrow=m, ncol=(k-3))
                            }
                            thita[,(k-2)]=runif(m, 0, dopi)
                            
                            fhi_gamma=output_fhi_s[t,]*vector_gamma(output_p_s[t,,], thita)
                            ############## etas:
                            fhi_eta=sqrt(1-output_fhi_s[t,]^2)*(vector_angles_nu(output_nu_s[t+1,,]))
                            ################# Hasting independent:
                            comp_mean=output_mu_s[t+1,]+(sweep(fhi_gamma,MARGIN=1,output_sigma_s[t+1,],`*`)/sqrt(output_p_s[t,,]))
                            comp_sigma=sweep(fhi_eta,MARGIN=1,output_sigma_s[t+1,],`*`)/sqrt(output_p_s[t,,])
                            
                            v=1
                            mix_prob=matrix(0,ncol=m, nrow=n)
                            repeat{
                                mixprob=apply(cbind(comp_mean[v,], comp_sigma[v,], output_p_s[t,,][v,]), 1, function(u)u[3]*dnorm(xobs, u[1], u[2]))
                                mixprob=rowSums(mixprob)
                                mix_prob[,v]=mixprob
                                if(v==m)
                                break
                                v=v+1
                            }
                            
                            propl=mix_prob
                            lpropl=log(propl)
                            slpropl=colSums(lpropl, na.rm=FALSE, dims=1)
                            proplik=slpropl
                            
                            log_ro_thit=betta*proplik-betta*curlik
                            Accept_thit=(log(runif(length(betta)))<(log_ro_thit))
                            
                            thi[Accept_thit,]=thita[Accept_thit,]
                            output_thita_s[t+1,,]=thi
                            
                            slcur[Accept_thit]=slpropl[Accept_thit]
                        }
                        ######################### Hasting independent proposal for fhi and tau:
                        degata2=output_fhi_s[t,]^2
                        curlik=slcur+dbeta(degata2, alpha1, alpha1, log=TRUE)
                        
                        alp=degata2*eps2+1
                        rgalp=(1-degata2)*eps2+1
                        prop3=rbeta(m, alp, rgalp)
                        A=sqrt(prop3)
                        if(k==2){
                            A=apply(cbind(-sqrt(prop3), sqrt(prop3)), 1, function(u)sample(u,1))
                        }
                        
                        fhi_gamma=A*vector_gamma(output_p_s[t,,], output_thita_s[t+1,,])
                        ############## etas:
                        fhi_eta=sqrt(1-prop3)*(vector_angles_nu(output_nu_s[t+1,,]))
                        
                        comp_mean=output_mu_s[t+1,]+(sweep(fhi_gamma,MARGIN=1,output_sigma_s[t+1,],`*`)/sqrt(output_p_s[t,,]))
                        comp_sigma=sweep(fhi_eta,MARGIN=1,output_sigma_s[t+1,],`*`)/sqrt(output_p_s[t,,])
                        
                        v=1
                        mix_prob=matrix(0,ncol=m, nrow=n)
                        repeat{
                            mixprob=apply(cbind(comp_mean[v,], comp_sigma[v,], output_p_s[t,,][v,]), 1, function(u)u[3]*dnorm(xobs, u[1], u[2]))
                            mixprob=rowSums(mixprob)
                            mix_prob[,v]=mixprob
                            if(v==m)
                            break
                            v=v+1
                        }
                        propl=mix_prob
                        lpropl=log(propl)
                        slpropl=colSums(lpropl, na.rm=FALSE, dims=1)
                        proplik=slpropl+dbeta(prop3, alpha1, alpha1, log=TRUE)
                        
                        alppdgt=prop3*eps2+1
                        alppdgtt=(1-prop3)*eps2+1
                        
                        log_ro_dir=betta*proplik-betta*curlik+dbeta(degata2, alppdgt, alppdgtt, log=TRUE)-dbeta(prop3, alp, rgalp, log=TRUE)
                        Accept_degata=(log(runif(length(betta)))<(log_ro_dir))
                        
                        fhi[Accept_degata]=A[Accept_degata]
                        output_fhi_s[t+1,]=fhi
                        
                        slcur[Accept_degata]=slpropl[Accept_degata]
                        post_dist_fhi[t,]=exp(slcur+dbeta(output_fhi_s[t+1,], alpha1, alpha1, log=TRUE))
                        naccept4[Accept_degata]=naccept4[Accept_degata]+1
                        na4[Accept_degata[1]]=na4[Accept_degata[1]]+1
                        
                        if(a>300){men_adapt_fhi=median(adapt_degata_ra[(a-10):a])
                            cond=(men_adapt_fhi>.4&men_adapt_fhi<.5)
                            adapt_condition_fhi[cond==TRUE]=FALSE
                        }
                        if(adapt_condition_fhi==TRUE){
                            ac_degata_ra=na4/50
                            adapt_degata_ra[a][a_adapt==50]=ac_degata_ra
                            ac_cond_degata1=(ac_degata_ra<.44)
                            ac_cond_degata2=(ac_degata_ra>.44)
                            
                            coc1=(eps2<=5)
                            coc2=(eps2>5)
                            eps2[ac_cond_degata1&coc2&a_adapt==50]=eps2[ac_cond_degata1&coc2&a_adapt==50]+5
                            eps2[ac_cond_degata2&coc2&a_adapt==50]=eps2[ac_cond_degata2&coc2&a_adapt==50]-5
                            eps2[ac_cond_degata1&coc1&a_adapt==50]=eps2[ac_cond_degata1&coc1&a_adapt==50]+5
                            eps2[ac_cond_degata2&coc1&a_adapt==50]=eps2[ac_cond_degata2&coc1&a_adapt==50]-eps2[ac_cond_degata2&coc1&a_adapt==50]/100
                            adapt_eps_degata[a][a_adapt==50]=eps2
                            a_fh[a_adapt==50]=a_fh[a_adapt==50]+1}else{adapt_degata_ra[a][a_adapt==50]=adapt_degata_ra[a-1][a_adapt==50]
                                adapt_eps_degata[a][a_adapt==50]=adapt_eps_degata[a-1][a_adapt==50]}
                            ################### Hasting independent proposal for p:
                            curlik=slcur+log(apply(output_p_s[t,,], 1, function(u)den_dirichlet(u, alpha0)))
                            
                            alp1=output_p_s[t,,]*eps1+1
                            prop4=matrix(0, nrow=m, ncol=k)
                            v=1
                            repeat{
                                prop4[v,]=rdirichlet(1, alp1[v,])
                                if(v==m)
                                break
                                v=v+1
                            }
                            alppt=alp1
                            alpp=prop4*eps1+1
                            log_dirp1=rep(0,m)
                            log_dirp2=rep(0,m)
                            for(l in 1:m){
                                log_dirp1[l]=log(apply(matrix(output_p_s[t,,][l,], nrow=1), 1, function(u)den_dirichlet(u, alpp[l,])))
                                log_dirp2[l]=log(apply(matrix(prop4[l,], nrow=1), 1, function(u)den_dirichlet(u, alppt[l,])))
                            }
                            
                            fhi_gamma=output_fhi_s[t+1,]*vector_gamma(prop4, output_thita_s[t+1,,])
                            ############## etas:
                            fhi_eta=sqrt(1-output_fhi_s[t+1,]^2)*(vector_angles_nu(output_nu_s[t+1,,]))
                            ###############
                            comp_mean=output_mu_s[t+1,]+(sweep(fhi_gamma,MARGIN=1,output_sigma_s[t+1,],`*`)/sqrt(prop4))
                            comp_sigma=sweep(fhi_eta,MARGIN=1,output_sigma_s[t+1,],`*`)/sqrt(prop4)
                            
                            v=1
                            mix_prob=matrix(0,ncol=m, nrow=n)
                            repeat{
                                mixprob=apply(cbind(comp_mean[v,], comp_sigma[v,], prop4[v,]), 1, function(u)u[3]*dnorm(xobs, u[1], u[2]))
                                mixprob=rowSums(mixprob)
                                mix_prob[,v]=mixprob
                                if(v==m)
                                break
                                v=v+1
                            }
                            
                            propl=mix_prob
                            lpropl=log(propl)
                            slpropl=colSums(lpropl, na.rm=FALSE, dims=1)
                            proplik=slpropl+log(apply(prop4, 1, function(u)den_dirichlet(u, alpha0)))
                            
                            log_ro_p=betta*proplik-betta*curlik+log_dirp1-log_dirp2
                            Accept_p=(log(runif(length(betta)))<(log_ro_p))
                            
                            pp[Accept_p,]=prop4[Accept_p,]
                            output_p_s[t+1,,]=pp
                            
                            slcur[Accept_p]=slpropl[Accept_p]
                            post_dist_p[t,]=exp(slcur+log(apply(output_p_s[t+1,,], 1, function(u)den_dirichlet(u, alpha0))))
                            naccept3[Accept_p]=naccept3[Accept_p]+1
                            na3[Accept_p[1]]=na3[Accept_p[1]]+1
                            
                            if(a>300){men_adapt_p=median(adapt_p_ra[(a-10):a])
                                if(k==2){cond=(men_adapt_p>.4&men_adapt_p<.5)}else{
                                    cond=(men_adapt_p>.2&men_adapt_p<.28)
                                }
                                adapt_condition_p[cond==TRUE]=FALSE
                            }
                            
                            if(adapt_condition_p==TRUE){
                                ac_p_ra=na3/50
                                adapt_p_ra[a][a_adapt==50]=ac_p_ra
                                
                                if(k==2){
                                    ac_cond_p1=(ac_p_ra<.44)
                                    ac_cond_p2=(ac_p_ra>.44)
                                    
                                    coco1=(eps1<=5)
                                    coco2=(eps1>5)
                                    eps1[ac_cond_p1&coco2&a_adapt==50]=eps1[ac_cond_p1&coco2&a_adapt==50]+5
                                    eps1[ac_cond_p2&coco2&a_adapt==50]=eps1[ac_cond_p2&coco2&a_adapt==50]-5
                                    eps1[ac_cond_p1&coco1&a_adapt==50]=eps1[ac_cond_p1&coco1&a_adapt==50]+5
                                    eps1[ac_cond_p2&coco1&a_adapt==50]=eps1[ac_cond_p2&coco1&a_adapt==50]-eps1[ac_cond_p2&coco1&a_adapt==50]/100}else{
                                        ac_cond_p1=(ac_p_ra<.234)
                                        ac_cond_p2=(ac_p_ra>.234)
                                        
                                        coco1=(eps1<=5)
                                        coco2=(eps1>5)
                                        eps1[ac_cond_p1&coco2&a_adapt==50]=eps1[ac_cond_p1&coco2&a_adapt==50]+5
                                        eps1[ac_cond_p2&coco2&a_adapt==50]=eps1[ac_cond_p2&coco2&a_adapt==50]-5
                                        eps1[ac_cond_p1&coco1&a_adapt==50]=eps1[ac_cond_p1&coco1&a_adapt==50]+5
                                        eps1[ac_cond_p2&coco1&a_adapt==50]=eps1[ac_cond_p2&coco1&a_adapt==50]-eps1[ac_cond_p2&coco1&a_adapt==50]/100}
                                    adapt_eps_p[a][a_adapt==50]=eps1
                                a_p[a_adapt==50]=a_p[a_adapt==50]+1}else{adapt_p_ra[a][a_adapt==50]=adapt_p_ra[a-1][a_adapt==50]
                                    adapt_eps_p[a][a_adapt==50]=adapt_eps_p[a-1][a_adapt==50]}
                                if(k>2){
                                    ######################### Hasting independent proposal for theta:
                                    curlik=slcur
                                    
                                    born1=output_thita_s[t+1,,]-eps
                                    born1=matrix(born1, nrow=m, ncol=(k-2))
                                    born2=output_thita_s[t+1,,]+eps
                                    born2=matrix(born2, nrow=m, ncol=(k-2))
                                    thita=matrix(0, nrow=m, ncol=(k-2))
                                    if(k>3){
                                        v=1
                                        repeat{
                                            tha=runif(m, born1[,v], born2[,v])
                                            co1=(tha<0)
                                            co2=(tha>pi)
                                            tha[co1]=pi+tha[co1]
                                            tha[co2]=tha[co2]-pi
                                            thita[,v]=tha
                                            if(v==(k-3))
                                            break
                                            v=v+1
                                        }
                                    }
                                    tha=runif(m, born1[,(k-2)], born2[,(k-2)])
                                    co1=(tha<0)
                                    co2=(tha>dopi)
                                    tha[co1]=dopi+tha[co1]
                                    tha[co2]=tha[co2]-dopi
                                    thita[,(k-2)]=tha
                                    fhi_gamma=output_fhi_s[t+1,]*vector_gamma(output_p_s[t+1,,], thita)
                                    ###############
                                    comp_mean=output_mu_s[t+1,]+(sweep(fhi_gamma,MARGIN=1,output_sigma_s[t+1,],`*`)/sqrt(output_p_s[t+1,,]))
                                    comp_sigma=sweep(fhi_eta,MARGIN=1,output_sigma_s[t+1,],`*`)/sqrt(output_p_s[t+1,,])
                                    
                                    v=1
                                    mix_prob=matrix(0,ncol=m, nrow=n)
                                    repeat{
                                        mixprob=apply(cbind(comp_mean[v,], comp_sigma[v,], output_p_s[t+1,,][v,]), 1, function(u)u[3]*dnorm(xobs, u[1], u[2]))
                                        mixprob=rowSums(mixprob)
                                        mix_prob[,v]=mixprob
                                        if(v==m)
                                        break
                                        v=v+1
                                    }
                                    
                                    propl=mix_prob
                                    lpropl=log(propl)
                                    slpropl=colSums(lpropl, na.rm=FALSE, dims=1)
                                    proplik=slpropl
                                    
                                    log_ro_thit=betta*proplik-betta*curlik
                                    Accept_thit=(log(runif(length(betta)))<(log_ro_thit))
                                    
                                    thi[Accept_thit,]=thita[Accept_thit,]
                                    output_thita_s[t+1,,]=thi
                                    
                                    slcur[Accept_thit]=slpropl[Accept_thit]
                                    post_dist_thita[t,]=exp(slcur-log(2*pi)-(k-3)*log(pi))
                                    naccept5[Accept_thit]=naccept5[Accept_thit]+1
                                    na5[Accept_thit[1]]=na5[Accept_thit[1]]+1
                                    
                                    if(a>300){men_adapt_thita=median(adapt_thita_ra[(a-10):a])
                                        if(k==3){cond=(men_adapt_thita>.4&men_adapt_thita<.5)}else{
                                            cond=(men_adapt_thita>.2&men_adapt_thita<.28)
                                        }
                                        adapt_condition_thita[cond==TRUE]=FALSE
                                    }
                                    
                                    if(adapt_condition_thita==TRUE){
                                        ac_thita_ra=na5/50
                                        adapt_thita_ra[a][a_adapt==50]=ac_thita_ra
                                        if(k==3){
                                            ac_cond_thita1=(ac_thita_ra<.44)
                                            ac_cond_thita2=(ac_thita_ra>.44)
                                            eps[ac_cond_thita1&a_adapt==50]=eps[ac_cond_thita1&a_adapt==50]-min(eps/100,eps/sqrt(t))
                                            eps[ac_cond_thita2&a_adapt==50]=eps[ac_cond_thita2&a_adapt==50]+min(eps/100,eps/sqrt(t))
                                        }else{
                                            ac_cond_thita1=(ac_thita_ra<.234)
                                            ac_cond_thita2=(ac_thita_ra>.234)
                                            eps[ac_cond_thita1&a_adapt==50]=eps[ac_cond_thita1&a_adapt==50]-min(eps/100,eps/sqrt(t))
                                            eps[ac_cond_thita2&a_adapt==50]=eps[ac_cond_thita2&a_adapt==50]+min(eps/100,eps/sqrt(t))
                                        }
                                        adapt_eps_thita[a][a_adapt==50]=eps
                                        a_theta[a_adapt==50]=a_theta[a_adapt==50]+1}else{adapt_thita_ra[a][a_adapt==50]=adapt_thita_ra[a-1][a_adapt==50]
                                            adapt_eps_thita[a][a_adapt==50]=adapt_eps_thita[a-1][a_adapt==50]}
                                }
                                
                                ######################### Hasting independent proposal for nus:
                                curlik=slcur
                                
                                born1=output_nu_s[t+1,,]-eps_nu
                                born1=matrix(born1, nrow=m, ncol=(k-1))
                                born2=output_nu_s[t+1,,]+eps_nu
                                born2=matrix(born2, nrow=m, ncol=(k-1))
                                nu_prop=matrix(0, nrow=m, ncol=(k-1))
                                sepido=3*pi/2
                                dopi=2*pi
                                if(k>2){
                                    v=1
                                    repeat{
                                        nha=runif(m, born1[,v], born2[,v])
                                        co1=(nha>pido & nha<=pi)
                                        co2=(nha>pi & nha<=sepido)
                                        co3=(nha>sepido & nha<=dopi)
                                        co4=(nha<0)
                                        nha[co1]=abs(nha[co1]-pi)
                                        nha[co2]=nha[co2]-pi
                                        nha[co3]=abs(nha[co3]-dopi)
                                        nha[co4]=abs(nha[co4])
                                        nu_prop[,v]=nha
                                        if(v==(k-2))
                                        break
                                        v=v+1
                                    }
                                }
                                nha=runif(m, born1[,(k-1)], born2[,(k-1)])
                                co1=(nha>pido & nha<=pi)
                                co2=(nha>pi & nha<=sepido)
                                co3=(nha>sepido & nha<=dopi)
                                co4=(nha<0)
                                nha[co1]=abs(nha[co1]-pi)
                                nha[co2]=nha[co2]-pi
                                nha[co3]=abs(nha[co3]-dopi)
                                nha[co4]=abs(nha[co4])
                                nu_prop[,(k-1)]=nha
                                
                                fhi_gamma=output_fhi_s[t+1,]*vector_gamma(output_p_s[t+1,,], output_thita_s[t+1,,])
                                ############## etas:
                                fhi_eta=sqrt(1-output_fhi_s[t+1,]^2)*(vector_angles_nu(nu_prop))
                                ################# Hasting independent:
                                comp_mean=output_mu_s[t+1,]+(sweep(fhi_gamma,MARGIN=1,output_sigma_s[t+1,],`*`)/sqrt(output_p_s[t+1,,]))
                                comp_sigma=sweep(fhi_eta,MARGIN=1,output_sigma_s[t+1,],`*`)/sqrt(output_p_s[t+1,,])
                                
                                if (sum(comp_sigma<0)==0){
                                    v=1
                                    mix_prob=matrix(0,ncol=m, nrow=n)
                                    repeat{
                                        mixprob=apply(cbind(comp_mean[v,], comp_sigma[v,], output_p_s[t+1,,][v,]), 1, function(u)u[3]*dnorm(xobs, u[1], u[2]))
                                        mixprob=rowSums(mixprob)
                                        mix_prob[,v]=mixprob
                                        if(v==m)
                                        break
                                        v=v+1
                                    }
                                    
                                    propl=mix_prob
                                    lpropl=log(propl)
                                    slpropl=colSums(lpropl, na.rm=FALSE, dims=1)
                                    proplik=slpropl
                                    
                                    log_ro_nu=betta*proplik-betta*curlik
                                    Accept_nu=(log(runif(length(betta)))<(log_ro_nu))
                                    
                                    nu[Accept_nu,]=nu_prop[Accept_nu,]
                                    output_nu_s[t+1,,]=nu
                                    
                                    slcur[Accept_nu]=slpropl[Accept_nu]
                                    post_dist_nu[t,]=exp(slcur-log(2*pi)-(k-2)*log(pi))
                                    naccept6[Accept_nu]=naccept6[Accept_nu]+1
                                    na6[Accept_nu[1]]=na6[Accept_nu[1]]+1
                                    
                                    if(a>300){men_adapt_nu=median(adapt_nu_ra[(a-10):a])
                                        if(k==2){cond=(men_adapt_nu>.4&men_adapt_nu<.5)}else{
                                            cond=(men_adapt_nu>.2&men_adapt_nu<.28)
                                        }
                                        adapt_condition_nu[cond==TRUE]=FALSE
                                    }
                                    
                                }else{
                                    output_nu_s[t+1,,]=output_nu_s[t,,]
                                }
                                
                                fhi_gamma=output_fhi_s[t+1,]*vector_gamma(output_p_s[t+1,,], output_thita_s[t+1,,])
                                ############## etas:
                                fhi_eta=sqrt(1-output_fhi_s[t+1,]^2)*(vector_angles_nu(output_nu_s[t+1,,]))
                                comp_mean=output_mu_s[t+1,]+(sweep(fhi_gamma,MARGIN=1,output_sigma_s[t+1,],`*`)/sqrt(output_p_s[t+1,,]))
                                comp_sigma=sweep(fhi_eta,MARGIN=1,output_sigma_s[t+1,],`*`)/sqrt(output_p_s[t+1,,])
                                component_mean[t+1,,]=comp_mean
                                component_sigma[t+1,,]=comp_sigma
                                
                                if(adapt_condition_nu==TRUE){
                                    ac_nu_ra=na6/50
                                    adapt_nu_ra[a][a_adapt==50]=ac_nu_ra
                                    if(k==2){
                                        ac_cond_nu1=(ac_nu_ra<.44)
                                        ac_cond_nu2=(ac_nu_ra>.44)
                                        eps_nu[ac_cond_nu1&a_adapt==50]=eps_nu[ac_cond_nu1&a_adapt==50]-min(eps_nu/100,eps_nu/sqrt(t))
                                        eps_nu[ac_cond_nu2&a_adapt==50]=eps_nu[ac_cond_nu2&a_adapt==50]+min(eps_nu/100,eps_nu/sqrt(t))
                                    }else{
                                        ac_cond_nu1=(ac_nu_ra<.234)
                                        ac_cond_nu2=(ac_nu_ra>.234)
                                        eps_nu[ac_cond_nu1&a_adapt==50]=eps_nu[ac_cond_nu1&a_adapt==50]-min(eps_nu/100,eps_nu/sqrt(t))
                                        eps_nu[ac_cond_nu2&a_adapt==50]=eps_nu[ac_cond_nu2&a_adapt==50]+min(eps_nu/100,eps_nu/sqrt(t))
                                    }
                                    adapt_eps_nu[a][a_adapt==50]=eps_nu
                                    a_nu[a_adapt==50]=a_nu[a_adapt==50]+1}else{adapt_nu_ra[a][a_adapt==50]=adapt_nu_ra[a-1][a_adapt==50]
                                        adapt_eps_nu[a][a_adapt==50]=adapt_eps_nu[a-1][a_adapt==50]}
                                    na1[a_adapt==50]=0
                                    na2[a_adapt==50]=0
                                    na3[a_adapt==50]=0
                                    na4[a_adapt==50]=0
                                    na5[a_adapt==50]=0
                                    na6[a_adapt==50]=0
                                    a[a_adapt==50]=a[a_adapt==50]+1
                                    a_adapt[a_adapt==50]=0
                                    stop_condition1=(adapt_condition_mu==FALSE & adapt_condition_sigma==FALSE & adapt_condition_fhi==FALSE & adapt_condition_p==FALSE & adapt_condition_thita==FALSE)
                                    if(rept==1 & stop_condition1==TRUE)
                                    break
                                    
                                    if(rept==2 & t==T_stop){
                                        criter1=rep(0,k)
                                        criter2=rep(0,k)
                                        T_stop1=T_stop-1000
                                        
                                        for(i in 1:k){
                                            criter1[i]=unlist(gelman.diag(mcmc.list(mcmc(component_mean[1:T_stop,1,i]), mcmc(component_mean[1:T_stop,2,i]))))[1]
                                            criter2[i]=unlist(gelman.diag(mcmc.list(mcmc(component_sigma[1:T_stop,1,i]), mcmc(component_sigma[1:T_stop,2,i]))))[1]
                                        }
                                        min_criter=min(c(criter1, criter2))
                                        max_criter=max(c(criter1, criter2))
                                        criter_condition=(min_criter>.9 & max_criter<3)
                                        T_stop=T_stop+1000
                                    }
                                    # stop_condition=(rept==2 & criter_condition==TRUE)
                                    #if(stop_condition==TRUE & t==(Nsim-1))
                                    #break
                                    stop_condition=(rept==2 & criter_condition==FALSE & t==(Nsim-1))
                                    
                                    if(stop_condition==TRUE){
                                        cat("WARNING! NOT CONVERGENT!", "\n")
                                    }
        }
        
        if(rept==2)
        break
        
        adapt_condition=FALSE
        adapt_condition_mu=FALSE
        adapt_condition_sigma=FALSE
        betta=rep(1, 4)
        T=Nsim
        m=length(betta)
        p0=output_p_s[t+1,,][1:m,]
        if(k>3){
            thita0=output_thita_s[t+1,,][1:m,]
        }else if(k==3){thita0=output_thita_s[t+1,,][1:m]}
        if(k>2){
            nu0=output_nu_s[t+1,,][1:m,]
        }else{
            nu0=output_nu_s[t+1,,][1:m]
        }
        fhi0=output_fhi_s[t+1,][1:m]
        naccept1=rep(0,m)
        naccept2=rep(0,m)
        naccept3=rep(0,m)
        naccept4=rep(0,m)
        naccept5=rep(0,m)
        naccept6=rep(0,m)
        
        n=length(xobs)
        meanobs=mean(xobs)
        
        mui=rep(mean(xobs), m)
        output_mu_s=matrix(mean(xobs), nrow=T, ncol=m)
        
        ss2=rep(var(xobs), m)
        output_sigma2_s=matrix(var(xobs), nrow=T, ncol=m)
        ss=rep(sd(xobs), m)
        output_sigma_s=sqrt(output_sigma2_s)
        
        pp=matrix(p0, nrow=m, ncol=k)
        output_p_s=array(p0, dim=c(T,m,k))
        for(i in 1:m){
            for(j in 1:k){
                output_p_s[,i,j]=p0[i,j]
            }
        }
        
        # if k>2
        if(k>3){
            thi=matrix(thita0, nrow=m, ncol=(k-2))
            output_thita_s=array(thita0, dim=c(T,m,(k-2)))
            for(i in 1:m){
                for(j in 1:(k-2)){
                    output_thita_s[,i,j]=thita0[i,j]
                }
            }
        }else if(k==3){
            thi=matrix(thita0, nrow=m, ncol=(k-2))
            output_thita_s=array(thita0, dim=c(T,m,(k-2)))
            for(i in 1:m){
                output_thita_s[,i,1]=thita0[i]
            }
        }else{output_thita_s=0}
        
        nu=matrix(nu0, nrow=m, ncol=(k-1))
        output_nu_s=array(nu0, dim=c(T,m,(k-1)))
        if(k>2){
            for(i in 1:m){
                for(j in 1:(k-1)){
                    output_nu_s[,i,j]=nu0[i,j]
                }
            }
        }else{
            for(i in 1:m){
                output_nu_s[,i,1]=nu0[i]
            }
        }
        
        fhi=fhi0
        output_fhi_s=matrix(fhi0, nrow=T, ncol=m, byrow=TRUE)
        component_mean=array(0, dim=c(T,m,k))
        component_sigma=array(0, dim=c(T,m,k))
        F=array(0, dim=c(m, k, (k-1)))
        post_dist_mu=matrix(0, nrow=T, ncol=m)
        post_dist_sigma=matrix(0, nrow=T, ncol=m)
        post_dist_p=matrix(0, nrow=T, ncol=m)
        post_dist_fhi=matrix(0, nrow=T, ncol=m)
        post_dist_thita=matrix(0, nrow=T, ncol=m)
        post_dist_nu=matrix(0, nrow=T, ncol=m)
        rept=rept+1
    }
    adapt_mu_ra=adapt_mu_ra[1:a_mu]
    adapt_sig_ra=adapt_sig_ra[1:a_sigma]
    adapt_p_ra=adapt_p_ra[1:a_p]
    adapt_degata_ra=adapt_degata_ra[1:a_fh]
    adapt_thita_ra=adapt_thita_ra[1:a_theta]
    adapt_nu_ra=adapt_nu_ra[1:a_nu]
    adapt_sd_mu=adapt_sd_mu[1:a_mu]
    adapt_sd_sig=adapt_sd_sig[1:a_sigma]
    adapt_eps_p=adapt_eps_p[1:a_p]
    adapt_eps_degata=adapt_eps_degata[1:a_fh]
    adapt_eps_thita=adapt_eps_thita[1:a_theta]
    adapt_eps_nu=adapt_eps_nu[1:a_nu]
    adapt_rat=list(adapt_mu_ra, adapt_sig_ra, adapt_p_ra, adapt_degata_ra, adapt_thita_ra, adapt_nu_ra)
    adapt_scale=list(adapt_sd_mu, adapt_sd_sig, adapt_eps_p, adapt_eps_degata, adapt_eps_thita, adapt_eps_nu)
    output_mu_s=output_mu_s[1:t,]
    output_sigma_s=output_sigma_s[1:t,]
    output_sigma2_s=output_sigma2_s[1:t,]
    output_fhi_s=abs(output_fhi_s[1:t,])
    output_p_s=output_p_s[1:t,,]
    output_nu_s=output_nu_s[1:t,,]
    component_mean=component_mean[1:t,,]
    component_sigma=component_sigma[1:t,,]
    n1=naccept1/t
    n2=naccept2/t
    n3=naccept3/t
    n4=naccept4/t
    n6=naccept6/t
    mins=rep(0,m)
    for(i in 1:m){
        mins[i]=min(n1[i], n2[i])
    }
    l=which.max(mins)
    mean_global=output_mu_s[,l] #[!output_mu_s[,l] %in% boxplot.stats(output_mu_s[,l])$out]
    sigma_global=output_sigma_s[,l] #[!output_sigma_s[,l] %in% boxplot.stats(output_sigma_s[,l])$out]
    phi=abs(output_fhi_s[,l])
    weights=output_p_s[,l,]
    epsilon_p=eps1
    epsilon_phi=eps2
    epsilon_xi=eps_nu
    epsilon_varpi=eps
    if(k==2){
        angles_xi=output_nu_s[1:t,l]
    }else{
        angles_xi=output_nu_s[1:t,l,]}
    component_means=component_mean[,l,]
    component_sigmas=component_sigma[,l,]
    
    if(k>2){
        angles_varpi=output_thita_s[1:t,l,]
        optimal_para=c(sd_mu, sd_sig, epsilon_p, epsilon_phi, epsilon_varpi, epsilon_xi)
        n5=naccept5/t
        accept_rat=c(n1[l], n2[l], n3[l], n4[l], n5[l], n6[l])
        estimate=list(mean_global, sigma_global, weights, angles_xi, phi, angles_varpi, accept_rat, optimal_para, adapt_rat, adapt_scale, component_means, component_sigmas)
        names(estimate)=c("mean_global", "sigma_global", "weights", "angles_xi", "phi", "angles_varpi", "accept_rat", "optimal_para", "adapt_rat", "adapt_scale", "component_means", "component_sigmas")
    }else{
        optimal_para=c(sd_mu, sd_sig, epsilon_p, epsilon_phi, epsilon_xi)
        accept_rat=c(n1[l], n2[l], n3[l], n4[l], n6[l])
        estimate=list(mean_global, sigma_global, weights, angles_xi, phi, accept_rat, optimal_para, adapt_rat, adapt_scale, component_means, component_sigmas)
        names(estimate)=c("mean_global", "sigma_global", "weights", "angles_xi", "phi", "accept_rat", "optimal_para", "adapt_rat", "adapt_scale", "component_means", "component_sigmas")
    }
    return(estimate)
}