init.p.trans <-
function(K,trans.const=TRUE)
{
    p.trans <- array(0,dim=c(K,K+1,K+1)) 
    if(K==1) p.trans[1,1,1] <- p.trans[1,1,2] <- p.trans[1,2,1] <- 1
    else
    {
        if(trans.const)
        {
            for(i in 1:K) for(j in 1:K) for (k in 1:K)
            {
                if(i!=j&&i!=k) p.trans[i,j,k] <- 0
                else
                {
                    if(i==j&&i==k)p.trans[i,j,k] <- 1
                    else p.trans[i,j,k] <- 0.5
                }
            }
            diag(p.trans[,,K+1]) <- 1
            diag(p.trans[,K+1,]) <- 1
        }
        else p.trans[1:K,1:(K+1),] <- p.trans[1:K,,1:(K+1)] <- 1/K
    }
    p.trans
}

