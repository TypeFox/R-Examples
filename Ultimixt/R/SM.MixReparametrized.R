SM.MixReparametrized <-function(xobs, estimate){
    k=dim(estimate[[3]])[2]
    mean_global=estimate[[1]][!estimate[[1]] %in% boxplot.stats(estimate[[1]])$out]
    sigma_global=estimate[[2]][!estimate[[2]] %in%boxplot.stats(estimate[[2]])$out]
    # Global mean:
    Mixture_mean=list(c(summary(mean_global)[3], summary(mean_global)[4]))
    Mixture_mean=data.frame(Mixture_mean)
    nameM=c("Mean: Mean of mixture distribution")
    names(Mixture_mean)=nameM
    # Global sigma:
    Mixture_sigma=list(c(summary(sigma_global)[3], summary(sigma_global)[4]))
    Mixture_sigma=data.frame(Mixture_sigma)
    nameS=c("Sd: Sd of mixture distribution")
    names(Mixture_sigma)=nameS
    # Phi:
    Mixture_phi=list(c(summary(estimate[[5]])[3], summary(estimate[[5]])[4]))
    Mixture_phi=data.frame(Mixture_phi)
    nameP=c("Phi")
    names(Mixture_phi)=nameP
    
    # theta:
    if(k>2){
        tta=cbind(c(estimate[[6]]))
        number=k-2
        results=kmeans(tta, number)
        angle1=list(c(results$centers))
        names(angle1)=c("Angles. 1.")
    }
    
    # xi:
    exi=cbind(c(estimate[[4]]))
    number=k-1
    results=kmeans(exi, number)
    angle2=list(c(results$centers))
    names(angle2)=c("Angles. 2.")
    
    # component means, standard deviations and weights:
    if(k>2){T=dim(estimate[[11]])[1]
        mus=estimate[[11]]
        sigmas=estimate[[12]]}else{T=dim(estimate[[10]])[1]
            mus=estimate[[10]]
            sigmas=estimate[[11]]
        }
        ps=estimate[[3]]
        
        perm_mus=matrix(0, nrow=T, ncol=k)
        perm_sigmas=matrix(0, nrow=T, ncol=k)
        perm_p=matrix(0, nrow=T, ncol=k)
        
        v=1
        t=1
        repeat{
            x=rbind(mus[t,], sigmas[t,], ps[t,])
            y=matrix(0, ncol=k, nrow=3)
            l=1
            repeat{
                MIN=which.min(x[3,])
                y[,l]=x[,MIN]
                x=x[,-MIN]
                cond=(is.vector(x)==TRUE)
                if(cond==TRUE)
                break
                l=l+1
            }
            y[,k]=x
            perm_mus[t,]=y[1,]
            perm_sigmas[t,]=y[2,]
            perm_p[t,]=y[3,]
            
            if(t==dim(mus)[1]-1)
            break
            t=t+1
        }
        
        mean_p=colSums(perm_p)/T
        f=1
        repeat{
            cond=(mean_p[f]>.06)
            if(cond==TRUE |f==k)
            break
            f=f+1
        }
        f=f-1
        
        k_prim=k-f
        new_perm_mus=perm_mus[1:t,(f+1):k]
        new_perm_sigmas=perm_sigmas[1:t,(f+1):k]
        
        dada=cbind(c(new_perm_mus), c(log(new_perm_sigmas)))
        results=kmeans(dada, k_prim)
        
        A=cbind(c(new_perm_mus), c(new_perm_sigmas), results$cluster, results$cluster)
        S=length(c(new_perm_mus))
        m1=s1=matrix(0, ncol=k_prim, nrow=S)
        means_cluster=sigmas_cluster=p_cluster=list(0)
        l2=rep(1,k_prim)
        l3=rep(1,k_prim)
        
        for(i in 1:S){
            j=1
            repeat{
                if(A[i,3]==j){m1[l2[j],j]=A[i,1];l2[j]=l2[j]+1}
                if(A[i,4]==j){s1[l3[j],j]=A[i,2];l3[j]=l3[j]+1}
                cond=(j==k_prim)
                if(cond==TRUE)
                break
                j=j+1
            }
        }
        l2=l2-1
        l3=l3-1
        for(j in 1:k_prim){
            means_cluster[[j]]=m1[1:l2[j],j]
            sigmas_cluster[[j]]=s1[1:l3[j],j]
        }
        
        summary_means_cluster=summary_sigmas_cluster=summary_p_cluster=list(0)
        if(f>0){
            for(j in 1:f){
                summary_means_cluster[[j]]=summary(perm_mus[1:t,j])[3:4]
                summary_sigmas_cluster[[j]]=summary(perm_sigmas[1:t,j])[3:4]
                summary_p_cluster[[j]]=summary(perm_p[1:t,j])[3:4]
            }
        }
        for(j in 1:k_prim){
            summary_means_cluster[[(f+j)]]=summary(means_cluster[[j]])[3:4]
            summary_sigmas_cluster[[(f+j)]]=summary(sigmas_cluster[[j]])[3:4]
            summary_p_cluster[[(f+j)]]=summary(perm_p[1:t,(f+j)])[3:4]
        }
        # Max likelihood:
        like_mean=rep(0, 10000)
        mean_weights=matrix(0, ncol=k, nrow=10000)
        mean_mus=matrix(0, ncol=k, nrow=10000)
        mean_sigs=matrix(0, ncol=k, nrow=10000)
        
        median_weights=matrix(0, ncol=k, nrow=10000)
        median_mus=matrix(0, ncol=k, nrow=10000)
        median_sigs=matrix(0, ncol=k, nrow=10000)
        x=seq(1,k,1)
        
        mean_weights[1,]=pps_mean=matrix(unlist(summary_p_cluster), ncol=2, byrow=TRUE)[,2]
        mean_mus[1,]=means_mean=matrix(unlist(summary_means_cluster), ncol=2, byrow=TRUE)[,2]
        mean_sigs[1,]=sigmas_mean=matrix(unlist(summary_sigmas_cluster), ncol=2, byrow=TRUE)[,2]
        median_weights[1,]=pps_median=matrix(unlist(summary_p_cluster), ncol=2, byrow=TRUE)[,1]
        median_mus[1,]=means_median=matrix(unlist(summary_means_cluster), ncol=2, byrow=TRUE)[,1]
        median_sigs[1,]=sigmas_median=matrix(unlist(summary_sigmas_cluster), ncol=2, byrow=TRUE)[,1]
        
        f_mean = rep(0,length(xobs),1)
        for (i in 1:length(pps_mean)){f_mean = f_mean + pps_mean[i]*dnorm(xobs,mean=means_mean[i],sd=sigmas_mean[i])}
        like_mean[1]=sum(f_mean)
        
        v=2
        repeat{
            f_mean = f_median=rep(0,length(xobs),1)
            y=permute(x)
            pps_mean=pps_mean[y]
            pps_median=pps_median[y]
            y=permute(x)
            means_mean=means_mean[y]
            means_median=means_median[y]
            y=permute(x)
            sigmas_mean=sigmas_mean[y]
            sigmas_median=sigmas_median[y]
            
            for (i in 1:length(pps_mean)){f_mean = f_mean + pps_mean[i]*dnorm(xobs,mean=means_mean[i],sd=sigmas_mean[i])}
            like_mean[v]=sum(f_mean)
            mean_weights[v,]=pps_mean
            mean_mus[v,]=means_mean
            mean_sigs[v,]=sigmas_mean
            
            median_weights[v,]=pps_median
            median_mus[v,]=means_median
            median_sigs[v,]=sigmas_median
            if(v==10000)
            break
            v=v+1
        }
        
        MT= which.max(like_mean)
        median_weights[MT,1]=1-sum(median_weights[MT,2:k])
        mean_weights[MT,1]=1-sum(mean_weights[MT,2:k])
        
        for(w in 1:k){
            summary_p_cluster[[w]]=c(median_weights[MT,w], mean_weights[MT,w])
            summary_means_cluster[[w]]=c(median_mus[MT,w], mean_mus[MT,w])
            summary_sigmas_cluster[[w]]=c(median_sigs[MT,w], mean_sigs[MT,w])
            
            names(summary_means_cluster[[w]])=c("Median", "mean")
            names(summary_sigmas_cluster[[w]])=c("Median", "mean")
            names(summary_p_cluster[[w]])=c("Median", "mean")
        }
        
        names(summary_means_cluster)=rep("mean",k)
        names(summary_sigmas_cluster)=rep("sd",k)
        names(summary_p_cluster)=rep("weight",k)
        
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
        
        print(Mixture_mean)
        cat("", iter <- c("##############################"), "\n")
        print(Mixture_sigma)
        cat("", iter <- c("##############################"), "\n")
        print(Mixture_phi)
        cat("", iter <- c("##############################"), "\n")
        if(k>2){
            print(angle1)
            cat("", iter <- c("##############################"), "\n")
        }
        print(angle2)
        cat("", iter <- c("##############################"), "\n")
        cat("Component means, standard deviations and weights:", iter <- c(""), "\n")
        cat("", iter <- c(""), "\n")
        print(data.frame(summary_p_cluster))
        cat("", iter <- c(""), "\n")
        print(data.frame(summary_means_cluster))
        cat("", iter <- c(""), "\n")
        print(data.frame(summary_sigmas_cluster))
        cat("", iter <- c("##############################"), "\n")
        print(Acc_rat)
        cat("", iter <- c("##############################"), "\n")
        print(Opt_scale)
}
