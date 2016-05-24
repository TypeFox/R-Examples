optim.noconst.ordi <-
function(y,status,weight,param,x=NULL,var.list=NULL)
{
    minp <- 0.00000001
    S <- apply(y,2,max,na.rm=TRUE)
    alpha <- list(NULL)
    for(j in 1:ncol(y))
    {
        if(!is.null(var.list[[j]])) S[j] <- S[j]+length(var.list[[j]])
        y[is.na(y[,j]),j] <- S[j]+10
        alpha[[j]] <- matrix(0,nrow=ncol(weight),ncol=S[j]-1)
        f.j.s.k <- t(apply(param$alpha[[j]],1,p.compute))
        for(k in 1:ncol(weight))
        {
            alpha[[j]][k,1] <- sum(weight[status==2,k]*(y[,j]==1))+sum(weight[status==0,k])*f.j.s.k[k,1]
            alpha[[j]][k,1] <- alpha[[j]][k,1]/sum(weight[,k])
            alpha[[j]][k,1] <- logit(min(max(alpha[[j]][k,1],minp),1-minp))
            if(S[j]>2) for(s in 2:(S[j]-1))
            {
                alpha[[j]][k,s] <- sum(weight[status==2,k]*(y[,j]<=s))+sum(weight[status==0,k])*sum(f.j.s.k[k,1:s])
                alpha[[j]][k,s] <- alpha[[j]][k,s]/sum(weight[,k])
                alpha[[j]][k,s] <- logit(min(max(alpha[[j]][k,s],minp),1-minp))-sum(alpha[[j]][k,1:(s-1)])
            }
        }
    }
    param <- list("alpha"=alpha)
    param
}

