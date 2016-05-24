optim.probs <-
function(ped,probs,optim.probs.indic=c(TRUE,TRUE,TRUE,TRUE),res.weight,famdep=TRUE)
{
    fam <- ped[,1]
    id <- ped[,2]
    dad <- ped[,3]
    mom <- ped[,4]
    status <- ped[,6]
    K <- length(probs$p)
    if(famdep)
    {
        found <- which(dad==0)
        children <- which(dad>0)
        connect <- intersect(children,union(which(id%in%dad),which(id%in%mom)))

        if(K==1) probs$p <- 1
        else
        {
            denom <- sum(res.weight$w[found,,1:K])
            if(denom>0) probs$p <- apply(res.weight$w[found,,1:K],3,sum)/denom
        }
        if(optim.probs.indic[1]) probs$p0 <- ifelse(sum(res.weight$w[found,2,])>0,sum(res.weight$w[found,2,K+1])/sum(res.weight$w[found,2,]),0)
        if(optim.probs.indic[2]&length(connect)>0)
        {
            mat.trio <- array(res.weight$ww[connect,2,,K+1,]+res.weight$ww[connect,2,,,K+1],dim=c(length(connect),K+1,K+1))
            probs$p0connect <- ifelse(apply(mat.trio,3,sum)>0,apply(matrix(mat.trio[,K+1,],nrow=length(connect),ncol=K+1),2,sum)/apply(mat.trio,3,sum),0)
            probs$p0connect <- probs$p0connect[1:K]
        }
        mat.trio <- array(res.weight$ww[children,1,,,]+res.weight$ww[children,2,,,],dim=c(length(children),rep(K+1,3)))
        mat.trio <- mat.trio+aperm(mat.trio,perm=c(1,2,4,3))
        for(c in 1:K)
        {
            denom <- apply(array(mat.trio[,1:K,,],dim=c(dim(mat.trio)[1],K,dim(mat.trio)[3],dim(mat.trio)[4])),c(3,4),sum)
            probs$p.trans[c,,] <- ifelse(denom>0,apply(mat.trio,c(2,3,4),sum)[c,,]/denom,0)
        }
        probs$p.trans[,K+1,K+1] <- 0
        if(optim.probs.indic[3]) probs$p.found <- sum(res.weight$w[found,1,1:K])/length(found)
        if(optim.probs.indic[4])
        {
            non.found <- (1:length(id))[-found]
            if(length(non.found)>0)
            {
                mat.trio <- array(res.weight$ww[non.found,1,1:K,,],dim=c(length(non.found),K,K+1,K+1))
                mat.trio <- mat.trio+aperm(mat.trio,perm=c(1,2,4,3))
                probs$p.child <- sum(mat.trio)/(2*length(non.found))
            }
        }
    }
    else
    {
        if(K>1)
        {
            denom <- sum(res.weight$w[,,1:K])
            if(denom>0) probs$p <- apply(res.weight$w[,,1:K],3,sum)/denom
        }
        else probs$p <- 1
        if(optim.probs.indic[1]) probs$p0 <- ifelse(sum(res.weight$w[,2,])>0,sum(res.weight$w[,2,K+1])/sum(res.weight$w[,2,]),0)
        if(optim.probs.indic[2]) probs$p.aff <- sum(res.weight$w[,1,1:K])/length(id)
    }
    probs
}

