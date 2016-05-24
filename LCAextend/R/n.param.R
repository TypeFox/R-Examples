n.param <-
function(y,K,trans.const=TRUE,optim.param,optim.probs.indic=c(TRUE,TRUE,TRUE,TRUE),famdep=TRUE)
{
    n.params <- K-1
    if(famdep)
    {
        if(optim.probs.indic[1]) n.param <- n.params+1
        if(optim.probs.indic[2]) n.param <- n.params+K
        if(trans.const) n.params <- n.params+K*(K-1)/2
        else n.params <- n.params+(K-1)*K*(K+3)/2
        if(optim.probs.indic[3]) n.param <- n.params+1
        if(optim.probs.indic[4]) n.param <- n.params+1
    }
    else if(optim.probs.indic[2]) n.param <- n.params+1

    if(identical(optim.param,optim.indep.norm)) n.params <- n.params+K*ncol(y)+K*ncol(y)
    if(identical(optim.param,optim.diff.norm)) n.params <- n.params+K*ncol(y)+ncol(y)*(ncol(y)+1)/2
    if(identical(optim.param,optim.equal.norm)) n.params <- n.params+K*ncol(y)+K*ncol(y)+K
    if(identical(optim.param,optim.gene.norm)) n.params <- n.params+K*ncol(y)+K*ncol(y)*(ncol(y)+1)/2
    n.cote <- NULL
    for(j in 1:ncol(y)) n.cote[j] <- length(table(y[,j]))
    if(identical(optim.param,optim.noconst.ordi)) n.params <- n.params+K*sum(n.cote-1)
    if(identical(optim.param,optim.const.ordi)) n.params <- n.params+ncol(y)*K+sum(n.cote-2)
    n.params
}

