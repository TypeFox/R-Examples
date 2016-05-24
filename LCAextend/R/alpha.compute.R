alpha.compute <-
function(p)
{
    alpha <- NULL
    if((!is.vector(p))|(any(p<0))|(any(p>1))|(abs(sum(p)-1)>.Machine$double.eps)) stop("p must be a vector of probabilities\n")

    if(length(p)==1) alpha <- NA
    else
    {
        if(p[1]==0) alpha[1] <- -Inf
        else if(p[1]==1) alpha[1] <- Inf
        else alpha[1] <- logit(p[1])
        for(i in setdiff(1:(length(p)-1),1))
        {
            cum.p <- sum(p[1:i])
            cum.alpha <- sum(alpha[!is.infinite(alpha)])
            alpha[i] <- ifelse(cum.p==1,Inf,logit(cum.p)-cum.alpha)
        }
    }
    alpha
  }

