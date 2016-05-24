p.post.child <-
function(child,c.connect,c.spouse,status,probs,fyc)
{
    p.child <- matrix(0,nrow=2,ncol=length(probs$p)+1)
    if(status[child]==2)
    {
        if((c.connect==length(probs$p)+1)&&(c.spouse==length(probs$p)+1)) p.child <- matrix(0,nrow=2,ncol=length(probs$p)+1)
        else
            {
                p.child[1,] <- fyc[child,]*c(probs$p.trans[,c.connect,c.spouse],0)*probs$p.child
                p.child[2,] <- 0
            }
    }
    if(status[child]==1)
    {
        if((c.connect<length(probs$p)+1)&&(c.spouse<length(probs$p)+1))
        {
            p.child[1,] <- 0
            p.child[2,] <- c(probs$p.trans[,c.connect,c.spouse],0)*(1-probs$p.child)
        }
        if((c.connect==length(probs$p)+1)&&(c.spouse==length(probs$p)+1))
        {
            p.child[1,] <- 0
            p.child[2,] <- c(rep(0,times=length(probs$p)+1-1),1)*(1-probs$p.child)
        }
        if((c.connect==length(probs$p)+1||c.spouse==length(probs$p)+1)&&min(c.connect,c.spouse)<length(probs$p)+1)
        {
            c.star <- min(c.connect,c.spouse)
            p.child[1,] <- 0
            p.child[2,] <- c((1-probs$p0connect[c.star])*probs$p.trans[,c.connect,c.spouse],probs$p0connect[c.star])*(1-probs$p.child)
        }
    }
    if(status[child]==0)
    {
        if((c.connect<length(probs$p)+1)&&(c.spouse<length(probs$p)+1))
        {
            p.child[1,] <- c(probs$p.trans[,c.connect,c.spouse],0)*probs$p.child
            p.child[2,] <- c(probs$p.trans[,c.connect,c.spouse],0)*(1-probs$p.child)
        }
        if((c.connect==length(probs$p)+1)&&(c.spouse==length(probs$p)+1))
        {
            p.child[1,] <- 0
            p.child[2,] <- c(rep(0,times=length(probs$p)+1-1),1)*(1-probs$p.child)
        }
        if((c.connect==length(probs$p)+1||c.spouse==length(probs$p)+1)&&min(c.connect,c.spouse)<length(probs$p)+1)
        {
            c.star <- min(c.connect,c.spouse)
            p.child[1,] <- c(probs$p.trans[,c.connect,c.spouse],0)*probs$p.child
            p.child[2,] <- c((1-probs$p0connect[c.star])*probs$p.trans[,c.connect,c.spouse],probs$p0connect[c.star])*(1-probs$p.child)
        }
    }
    p.child
}

