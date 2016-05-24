SM.MAP.MixReparametrized <-
function(estimate,xobs,alpha0,alpha){
    k=dim(estimate[[3]])[2]
    mean_global=estimate[[1]][!estimate[[1]] %in% boxplot.stats(estimate[[1]])$out]
    sigma_global=estimate[[2]][!estimate[[2]] %in%boxplot.stats(estimate[[2]])$out]
    # Global mean:
    Mixture_mean=(c(mean(mean_global), median(mean_global),quantile(mean_global,prob=c(0.025,0.975))))
    names(Mixture_mean)[1:2]<-c("Mean","Median")
    # Global sigma:
    Mixture_sigma=(c(mean(sigma_global), median(sigma_global),quantile(sigma_global,prob=c(0.025,0.975))))
    names(Mixture_sigma)[1:2]<-c("Mean","Median")
    # Phi:
    Mixture_phi=(c(mean(estimate[[5]]), median(estimate[[5]]),quantile(estimate[[5]],prob=c(0.025,0.975))))
    names(Mixture_phi)[1:2]<-c("Mean","Median")
    
    # component means, standard deviations and weights:
    if(k>2){T=dim(estimate[[11]])[1]
        mus=estimate[[11]]
        sigmas=estimate[[12]]}else{T=dim(estimate[[10]])[1]
            mus=estimate[[10]]
            sigmas=estimate[[11]]
        }
        ps=estimate[[3]]
        
        ### log-likelihood + log-prior
        ### Note that uniform priors for angular components are discarded.
        logL_mx=matrix(0,dim(mus)[1],length(xobs))
        for (i in 1:length(xobs)){ for(j in 1:k){logL_mx[,i]=logL_mx[,i]+ps[,j]*dnorm(xobs[i],mus[,j],sigmas[,j])}}
        logPost=rowSums(log(logL_mx))-log(estimate[[2]])+rowSums(alpha0*log(ps))+lgamma(k*alpha0)-k*lgamma(alpha0)+dbeta(estimate$phi^2,alpha,alpha,log=T);
        max_ind=which(logPost==max(logPost)); max_ind=max_ind[1];
        
        vertex_mu=t(replicate(T,mus[max_ind,])); vertex_sigma=t(replicate(T,sigmas[max_ind,]));
        vertex_p=t(replicate(T,ps[max_ind,]));
        ind_perm=permutations(k,k); dist=matrix(,T,dim(ind_perm)[1])
        
        for(i in 1:dim(ind_perm)[1]){
            dist[,i]=sqrt(rowSums((mus[,ind_perm[i,]]-vertex_mu)^2+(sigmas[,ind_perm[i,]]-vertex_sigma)^2+(ps[,ind_perm[i,]]-vertex_p)^2))}
        
        perm_mus= perm_sigmas = perm_p = matrix(0, nrow=T, ncol=k)
        perm_xi=matrix(0, nrow=T, ncol=k-1)
        if (k>2){perm_varpi=matrix(0, nrow=T, ncol=k-2)}
        
        for (i in 1:T){final_ind=which(dist[i,]==min(dist[i,]))
            perm_mus[i,]=mus[i,ind_perm[final_ind,]]
            perm_sigmas[i,]=sigmas[i,ind_perm[final_ind,]]
            perm_p[i,]=ps[i,ind_perm[final_ind,]]}
        
        summary_means_cluster=summary_sigmas_cluster=summary_p_cluster=list(0)
        
        for(j in 1:k){
            summary_means_cluster[[j]]= c(mean(perm_mus[,j]),median(perm_mus[,j]),quantile(perm_mus[,j],prob=c(0.025,0.975)))
            names(summary_means_cluster[[j]])[1:2]<-c("Mean","Median")
            summary_sigmas_cluster[[(j)]]=c(mean(perm_sigmas[,j]),median(perm_sigmas[,j]),quantile(perm_sigmas[,j],prob=c(0.025,0.975)))
            names(summary_sigmas_cluster[[j]])[1:2]<-c("Mean","Median")
            summary_p_cluster[[(j)]]=c(mean(perm_p[,j]),median(perm_p[,j]),quantile(perm_p[,j],prob=c(0.025,0.975)))
            names(summary_p_cluster[[j]])[1:2]<-c("Mean","Median")
        }
        
        names(summary_means_cluster)=rep("mean",k)
        names(summary_sigmas_cluster)=rep("sd",k)
        names(summary_p_cluster)=rep("weight",k)
        
        
        ####################### inverse eta
        find_xi <- function(sig,phi,global_sig,p){
            k=length(sig); xi=rep(0,k-1);
            xi[1]=acos(sig[1]/sqrt(1-phi^2)/global_sig*sqrt(p[1])); cumSin=sin(xi[1]);
            if(k>2){for (i in 2:(k-1)){xi[i]=acos(sig[i]/sqrt(1-phi^2)/global_sig/cumSin*sqrt(p[i])); cumSin=cumSin*sin(xi[i])}}
            return(xi) }
        
        ####################### inverse gam
        find_varpi <- function(mu,phi,global_sig,global_mu,p){
            k=length(mu); varpi=rep(0,k-2);
            v=matrix(0,(k-1),k); #orthonormal vector
            v[1,1]=-sqrt(p[2]); v[1,2]=sqrt(p[1])
            if (k>2){for (i in 2:(k-1)){
                for (j in 1:i){ v[i,j]=-sqrt(p[j]*p[(i+1)])/sqrt(sum(p[1:i]))}
                v[i,(j+1)]=sqrt(sum(p[1:i])) }}
            for (i in 1:(k-1)){ v[i,]=v[i,]/(replicate(k,sqrt(sum(v[i,]^2)))) }
            z=(mu-global_mu)/global_sig*sqrt(p)/phi
            ## c=coefficients of v, (k-1)-orthonormal vectors
            c=rep(0,(k-1)); c[(k-1)]=z[k]/v[(k-1),k]
            for (i in 2:(k-1)){
                c[(k-i)]=(z[(k-i+1)]-c[(k-i+1):(k-1)]%*%v[(k-i+1):(k-1),(k-i+1)])/v[(k-i),(k-i+1)]
            }
            varpi=rep(0,k-2)
            varpi[1]=acos(c[1]); cumSin=sin(varpi[1]);
            if (k>3){ for (i in 2:(k-2)){varpi[i]=acos(c[i]/cumSin); cumSin=cumSin*sin(varpi[i])}}
            return(varpi) }
        
        
        #### computing angle components corresponding to perm_mus, perm_sigmas, and perm_p
        
        if (k==3){for (i in 1:T){perm_xi[i,] <- find_xi(perm_sigmas[i,],estimate$phi[i],estimate$sigma_global[i],perm_p[i,])
            perm_varpi[i] = find_varpi(perm_mus[i,],estimate$phi[i],estimate$sigma_global[i],estimate$mean_global[i], perm_p[i,])}
        }else if(k==2){
            for (i in 1:T){perm_xi[i] <- find_xi(perm_sigmas[i,],estimate$phi[i],estimate$sigma_global[i],perm_p[i,])}
        }else{
            for (i in 1:T){perm_xi[i,] <- find_xi(perm_sigmas[i,],estimate$phi[i],estimate$sigma_global[i],perm_p[i,])
                perm_varpi[i,] = find_varpi(perm_mus[i,],estimate$phi[i],estimate$sigma_global[i],estimate$mean_global[i], perm_p[i,])}}
        
        summary_ang_sig_cluster=summary_ang_mu_cluster=list(0)
        perm_xi=as.matrix(perm_xi)
        for(j in 1:(k-1)){
            summary_ang_sig_cluster[[j]]= c(mean(perm_xi[,j]),median(perm_xi[,j]),quantile(perm_xi[,j],prob=c(0.025,0.975)))
            names(summary_ang_sig_cluster[[j]])[1:2]<-c("Mean","Median")
        }
        names(summary_ang_sig_cluster)=rep("angle_sigma",(k-1))
        
        if (k>2){perm_varpi=as.matrix(perm_varpi);
            for (j in 1:(k-2)){ summary_ang_mu_cluster[[(j)]]=c(mean(perm_varpi[,j]),median(perm_varpi[,j]),quantile(perm_varpi[,j],prob=c(0.025,0.975)))
                names(summary_ang_mu_cluster[[j]])[1:2]<-c("Mean","Median")
            }
            names(summary_ang_mu_cluster)=rep("angle_mu",(k-2))}
        
        
        ####determining the class of labels and count a number of label switching occurance
        d_perm=matrix(,T,factorial(k)); class=rep(0,T)
        MAPmu=mus[max_ind,]; MAPsigma=sigmas[max_ind,]; MAPps=ps[max_ind,];
        
        for(i in 1:dim(ind_perm)[1]){
            d_perm[,i]=sqrt(rowSums((mus-t(replicate(T, MAPmu[ind_perm[i,]])))^2+(sigmas-t(replicate(T,MAPsigma[ind_perm[i,]])))^2+(ps-t(replicate(T,MAPps[ind_perm[i,]])))^2))}
        class[1]=which(d_perm[1,]==min(d_perm[1,]));
        change=l_stay=0; stay=1;
        for (i in 2:T){ class[i]=which(d_perm[i,]==min(d_perm[i,])); if (class[i]==class[i-1]){stay=stay+1}else{change=change+1;l_stay=c(l_stay,stay); stay=1;}}
        l_stay=l_stay[-1]
        
        # Acceptance rate of proposals and optimal scales:
        if(k>2){
            Acc_rat=estimate[[7]]
            nameA=c("mu", "sigma", "p", "phi", "theta", "xi")
            names(Acc_rat)=nameA
            Acc_rat=list(Acc_rat)
            nameA=c("Acceptance rate of proposals")
            names(Acc_rat)=nameA
            
            Opt_scale=estimate[[8]]
            nameO=c("s_mu", "s_sigma", "epsilon_p", "epsilon_phi", "epsilon_theta", "epsilon_xi")
            names(Opt_scale)=nameO
            Opt_scale=list(Opt_scale)
            nameO=c("Optimal proposal scales")
            names(Opt_scale)=nameO
        }else{
            Acc_rat=estimate[[6]]
            nameA=c("mu", "sigma", "p", "phi", "xi")
            names(Acc_rat)=nameA
            Acc_rat=list(Acc_rat)
            nameA=c("Acceptance rate of proposals")
            names(Acc_rat)=nameA
            
            Opt_scale=estimate[[7]]
            nameO=c("s_mu", "s_sigma", "epsilon_p", "epsilon_phi", "epsilon_xi")
            names(Opt_scale)=nameO
            Opt_scale=list(Opt_scale)
            nameO=c("Optimal proposal scales")
            names(Opt_scale)=nameO
        }
        
        cat("", iter <- c("##############################"), "\n")
        cat("Global mean ", iter <- c(""), "\n")
        print(Mixture_mean)
        cat("", iter <- c("##############################"), "\n")
        cat("Global standard deviation ", iter <- c(""), "\n")
        print(Mixture_sigma)
        cat("", iter <- c("##############################"), "\n")
        cat("Radius(phi)", iter <- c(""), "\n")
        print(Mixture_phi)
        #cat("", iter <- c("##############################"), "\n")
        #if(k>2){
        #cat("Angular component (w)", iter <- c(""), "\n")
        #print(angle1)
        #cat("", iter <- c("##############################"), "\n")
        #}
        #cat("Angular component (eta)", iter <- c(""), "\n")
        #print(angle2)
        cat("", iter <- c("##############################"), "\n")
        cat("Component means, standard deviations and weights:", iter <- c(""), "\n")
        cat("", iter <- c(""), "\n")
        print(data.frame(summary_p_cluster))
        cat("", iter <- c(""), "\n")
        print(data.frame(summary_means_cluster))
        cat("", iter <- c(""), "\n")
        print(data.frame(summary_sigmas_cluster))
        cat("", iter <- c("##############################"), "\n")
        cat("Angle components associated with component means and standard deviations:", iter <- c(""), "\n")
        cat("", iter <- c(""), "\n")
        print(data.frame(summary_ang_sig_cluster))
        if (k>2){
            cat("", iter <- c(""), "\n")
            print(data.frame(summary_ang_mu_cluster))
        }
        
        print(Acc_rat)
        cat("", iter <- c("##############################"), "\n")
        print(Opt_scale) 
        
        if (k>2){
            result=list(perm_mus,perm_sigmas,perm_p,perm_xi,perm_varpi,Mixture_mean, Mixture_sigma, Mixture_phi, summary_means_cluster, summary_sigmas_cluster, summary_p_cluster,l_stay,change)
            names(result)=c("MU","SIGMA","P","Ang_SIGMA","Ang_MU","Global_Mean","Global_Std","Phi","component_mu","component_sigma","component_p","l_stay","n_switch")}else{result=list(perm_mus,perm_sigmas,perm_p,perm_xi,Mixture_mean, Mixture_sigma, Mixture_phi, summary_means_cluster, summary_sigmas_cluster, summary_p_cluster,l_stay,change)
                names(result)=c("MU","SIGMA","P","Ang_SIGMA","Global_Mean","Global_Std","Phi","component_mu","component_sigma","component_p","l_stay","n_switch")
                
            }
            
            return(result)
}
