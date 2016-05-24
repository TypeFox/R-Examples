###########################################################
##### Fonction pour calculer la probabilite a posteriori. #
##########################################################

### Prototype : post(X)

post=function(X,...)
{
    UseMethod("post")

}

### Prototype : entropy(X)
entropy=function(X,...)
{
    UseMethod("entropy")
}



post.mmeln=function(X,...,mu=X$param$mu,tau=X$param$tau,sigma=X$param$sigma)
{
    post=matrix()
    if(is.null(tau))
    {
        post=matrix(1,X$N,X$G)
    }
    else
    {
        eta=X$Z%*%matrix(tau,ncol=X$G-1)
        P=cbind(1,exp(eta))/apply(cbind(1,exp(eta)),1,sum)
        fd=numeric()
        for(i in 1:X$G)
        {
            if(!X$equalcov)
                Cov=cov.tsf(sigma[[i]],X$cov,X$p)
            else
                Cov=cov.tsf(sigma[[1]],X$cov,X$p)
            pred=X$Xg[[i]]%*%mu[[i]]
            fd=cbind(fd,P[,i]*dmnorm(X$Y,as.vector(pred),Cov))
        }
        post=fd/apply(fd,1,sum)
    }
    post
}


entropy.mmeln=function(X,...)
{
    tauij=post(X)
    -sum(tauij*log(tauij))
}

