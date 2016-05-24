p.post.found <-
function(found,status,probs,fyc)
{
    p.found <- matrix(0,nrow=2,ncol=length(probs$p)+1)
    if(status[found]==2)
    {
        p.found[1,] <- fyc[found,]*c(probs$p,0)*probs$p.found
        p.found[2,] <- 0
    }
    if(status[found]==1)
    {
        p.found[1,] <- 0
        p.found[2,] <- c((1-probs$p0)*probs$p,probs$p0)*(1-probs$p.found)
    }
    if(status[found]==0)
    {
        p.found[1,] <- c(probs$p,0)*probs$p.found
        p.found[2,] <- c((1-probs$p0)*probs$p,probs$p0)*(1-probs$p.found)
    }
    p.found
}

