weight.famdep <-
function(id,dad,mom,status,probs,fyc,peel)
{
    id.origin <- id
    standard <- function(vec) ifelse(vec%in%id.origin,which(id.origin==vec),0)
    id <- apply(t(id),2,standard)
    dad <- apply(t(dad),2,standard)
    mom <- apply(t(mom),2,standard)
    peel$couple <- cbind(apply(t(peel$couple[,1]),2,standard),apply(t(peel$couple[,2]),2,standard))
    for(generat in 1:peel$generation) peel$peel.connect[generat,] <- apply(t(peel$peel.connect[generat,]),2,standard)

    res.upward <- upward(id,dad,mom,status,probs,fyc,peel)
    res.downward <- downward(id,dad,mom,status,probs,fyc,peel,res.upward)
    
    #here we permute ww to standardize it like ww[child,s_child,c_child,s_parent1,c_parent1,s_parent2,c_parent2]
    ww <- aperm(res.downward$ww,perm=c(1,4,7,2,5,3,6))
    #sum over s_parent1 and s_parent 2 to get ww[i,s_i,c_i,c_1,c_2]
    ww <- apply(ww,c(1,2,3,5,7),sum)
    #normalizing ww and w by likelihood
    ll <- sum(ww)/sum(dad>0)
    ww <- ww/ll
    w <- res.downward$w/ll
    res <- list("ww"=ww,"w"=w,"ll"=log(ll))
    res
}

