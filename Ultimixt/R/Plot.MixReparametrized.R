Plot.MixReparametrized <-
function(xobs, estimate){
    k=dim(estimate[[3]])[2]
    n=length(xobs)
    # Estimated mixture density:
    if(k>2){s=dim(estimate[[11]])[1]
        ss=s-10000
        ss[ss<0]=1
        sm=s-ss
        mu=estimate[[11]][ss:s,]
        sig=estimate[[12]][ss:s,]}else{s=dim(estimate[[10]])[1]
            ss=s-10000
            ss[ss<0]=1
            sm=s-ss
            mu=estimate[[10]][ss:s,]
            sig=estimate[[11]][ss:s,]
        }
        p=estimate[[3]][ss:s,]
        mix=list(k=k, p=p, mu=mu, sig=sig)
        yt=list(0)
        x=seq(min(xobs),max(xobs),length=1000)
        for(j in 1:sm){
            y=mix$p[j,1]*dnorm(x,mean=mix$mu[j,1],sd=sqrt(mix$sig[j,1]))
            for (i in 2:mix$k) y=y+mix$p[j,i]*dnorm(x,mean=mix$mu[j,i],sd=sqrt(mix$sig[j,i]))
            yt[[j]]=y
        }
        mean_density=rep(0, 1000)
        men=rep(0,sm)
        for(i in 1:1000){
            for(j in 1:sm){
                men[j]=yt[[j]][i]
            }
            mean_density[i]=mean(men)
        }
        hist(xobs,prob=T,main="Fitted mixture density",xlab="",col="gold", nclass=100)
        if(length(yt)>500){
            for(i in (length(yt)-500):length(yt)){lines(x, yt[[i]], lty=2, col="chocolate")}
        }else{
            for(i in 1:length(yt)){lines(x, yt[[i]], lty=2, col="chocolate")}
        }
        lines(x, mean_density, lwd=2)
        
        # mean, sd, phi:
        dev.new();par(mfrow=c(1,3))
        boxplot(estimate[[1]], col="gold2", main=expression(mu), outline=FALSE)
        boxplot(estimate[[2]], col="gold2", main=expression(sigma), outline=FALSE)
        boxplot(estimate[[5]], col="gold2", main=expression(varphi), outline=FALSE)
        
        # evolustion of proposal scales and acceptarnce rates:
        if(k>2){
            dev.new(); par(mfcol=c(3,4))
            plot(estimate[[9]][[1]], type="l", ylab=expression(mu), main="adaptive acceptance rat")
            plot(estimate[[9]][[2]], type="l", ylab=expression(sigma), main="")
            plot(estimate[[9]][[3]], type="l", ylab=expression(p), main="")
            plot(estimate[[9]][[4]], type="l", ylab=expression(phi), main="adaptive acceptance rat")
            plot(estimate[[9]][[5]], type="l", ylab=expression(varpi), main="")
            plot(estimate[[9]][[6]], type="l", ylab=expression(xi), main="")
            plot(estimate[[10]][[1]], type="l", ylab=expression(mu), main=expression(s[mu]))
            plot(estimate[[10]][[2]], type="l", ylab=expression(sigma), main=expression(s[sigma]))
            plot(estimate[[10]][[3]], type="l", ylab=expression(p), main=expression(epsilon[p]))
            plot(estimate[[10]][[4]], type="l", ylab=expression(phi), main=expression(epsilon[varphi]))
            plot(estimate[[10]][[5]], type="l", ylab=expression(varpi), main=expression(epsilon[omega]))
            plot(estimate[[10]][[6]], type="l", ylab=expression(xi), main=expression(epsilon[xi]))
        }else{
            dev.new(); par(mfcol=c(2,5))
            plot(estimate[[8]][[1]], type="l", ylab=expression(mu), main="adaptive acceptance rat")
            plot(estimate[[8]][[2]], type="l", ylab=expression(sigma), main="")
            plot(estimate[[8]][[3]], type="l", ylab=expression(p), main="")
            plot(estimate[[8]][[4]], type="l", ylab=expression(phi), main="adaptive acceptance rat")
            plot(estimate[[8]][[6]], type="l", ylab=expression(xi), main="")
            plot(estimate[[9]][[1]], type="l", ylab=expression(mu), main=expression(s[mu]))
            plot(estimate[[9]][[2]], type="l", ylab=expression(sigma), main=expression(s[sigma]))
            plot(estimate[[9]][[3]], type="l", ylab=expression(p), main=expression(epsilon[p]))
            plot(estimate[[9]][[4]], type="l", ylab=expression(phi), main=expression(epsilon[varphi]))
            plot(estimate[[9]][[6]], type="l", ylab=expression(xi), main=expression(epsilon[xi]))
        }
        
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
            
            dev.new()
            if(k==2){par(mfrow=c(1,2))}else if(k==3){par(mfrow=c(1,3))}else if(k==4){par(mfrow=c(2,2))}else if(k==5){par(mfrow=c(1,5))}else if(k==6){par(mfrow=c(2,3))}else if(k==7){par(mfrow=c(1,7))}else if(k==8){par(mfrow=c(2,4))}else if(k==9| k==10){par(mfrow=c(2,5))}else if(k==11| k==12){par(mfrow=c(2,6))}else if(k==13| k==14){par(mfrow=c(2,7))}else if(k==15| k==16){par(mfrow=c(2,8))}
            
            if(f>0){
                for(i in 1:f){
                    boxplot(perm_mus[1:t,i], col="gold4", main=expression(mu[i]))
                }
            }
            
            for(j in 1:k_prim){
                l=f+j
                boxplot(means_cluster[[j]], outline=FALSE, col="gold4", main=expression(mu[l]))
            }
            
            dev.new()
            if(k==2){par(mfrow=c(1,2))}else if(k==3){par(mfrow=c(1,3))}else if(k==4){par(mfrow=c(2,2))}else if(k==5){par(mfrow=c(1,5))}else if(k==6){par(mfrow=c(2,3))}else if(k==7){par(mfrow=c(1,7))}else if(k==8){par(mfrow=c(2,4))}else if(k==9| k==10){par(mfrow=c(2,5))}else if(k==11| k==12){par(mfrow=c(2,6))}else if(k==13| k==14){par(mfrow=c(2,7))}else if(k==15| k==16){par(mfrow=c(2,8))} 
            
            if(f>0){
                for(i in 1:f){
                    boxplot(perm_sigmas[1:t,i], col="gold4", main=expression(sigma[i]))
                }
            }
            
            for(j in 1:k_prim){
                l=f+j
                boxplot(sigmas_cluster[[j]], outline=FALSE, col="gold4", main=expression(sigma[l]))
            } 
            
            dev.new()
            if(k==2){par(mfrow=c(1,2))}else if(k==3){par(mfrow=c(1,3))}else if(k==4){par(mfrow=c(2,2))}else if(k==5){par(mfrow=c(1,5))}else if(k==6){par(mfrow=c(2,3))}else if(k==7){par(mfrow=c(1,7))}else if(k==8){par(mfrow=c(2,4))}else if(k==9| k==10){par(mfrow=c(2,5))}else if(k==11| k==12){par(mfrow=c(2,6))}else if(k==13| k==14){par(mfrow=c(2,7))}else if(k==15| k==16){par(mfrow=c(2,8))} 
            
            for(j in 1:k){
                boxplot(perm_p[1:t,j], outline=FALSE, col="gold4", main=expression(p[j]))
            }
            
}
