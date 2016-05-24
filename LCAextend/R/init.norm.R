init.norm <-
function(y,K,x=NULL,var.list=NULL)
{
    mu <- sigma <- list()
    if(nrow(y)<K) stop("there are a few individuals to perform a model with ",K," classes\n")
    if(K==1)
    {
        mu[[1]] <- apply(y,2,mean,na.rm=TRUE)
        sigma[[1]] <- var(y,na.rm=TRUE)
    }
    else
    {
        y.complete <- y[!apply(t(apply(y,1,is.na)),1,all),]
        if(ncol(y)==1) y.complete <- matrix(y.complete,nrow=length(y),ncol=1)
        clust <- cutree(hclust(dist(y.complete,"euclidean")),K)
        for(i in which(table(clust)>1))
        {
            y.complete.clusti <- matrix(y.complete[clust==i,],nrow=sum(clust==i),ncol=ncol(y))
            mu[[i]] <- apply(y.complete.clusti,2,mean,na.rm=TRUE)
            sigma[[i]] <- var(y.complete.clusti,na.rm=TRUE)
            if(det(sigma[[i]])<0.001) sigma[[i]] <- var(y.complete,na.rm=TRUE)
        }
        for(i in which(table(clust)==1))
        {
            mu[[i]] <- y.complete[clust==i,]
            sigma[[i]] <- sigma[[which(sum(clust)>1)[1]]]
        }
        order.mu <- order(unlist(lapply(mu,mean)))
        mu <- mu[order.mu]
        sigma <- sigma[order.mu]
    }
    res <- list("mu"=mu,"sigma"=sigma)
    res
}

