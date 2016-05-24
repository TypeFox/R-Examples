e.step <-
function(ped,probs,param,dens,peel,x=NULL,var.list=NULL,famdep=TRUE)
{
    fam <- ped[,1]
    id <- ped[,2]
    dad <- ped[,3]
    mom <- ped[,4]
    sex <- ped[,5]
    status <- ped[,6]

	y <- as.matrix(ped[,-(1:6)])

    K <- length(probs$p)

    fyc <- matrix(1,nrow=length(id),ncol=K+1)
    if(any(status==2))
    {
		y.x.aff <- as.matrix(y[status==2,])
		if(!is.null(x)) y.x.aff <- cbind(y.x.aff,as.matrix(x[status==2,]))
        fyc[status==2,1:K] <- t(apply(y.x.aff,1,dens,param,var.list))

		if(K>1) probi0 <- apply(as.matrix(fyc[status==2,1:K]),1,function(vec) all(vec<.Machine$double.eps))
        else probi0 <- fyc[status==2,1]<.Machine$double.eps
        if(any(probi0))
        {
            fyc[status==2,1:K][probi0] <- .Machine$double.eps
            warning("density less than ",round(.Machine$double.eps,18)," for at least one observation\n")
        }
    }
    w <- array(0,dim=c(length(id),2,K+1))
    ll <- 0
    if(famdep)
    {
        ww <- array(NA,dim=c(length(id),2,K+1,K+1,K+1))
        for(f in unique(fam))
        {
            res <- weight.famdep(id[fam==f],dad[fam==f],mom[fam==f],status[fam==f],probs,fyc[fam==f,],peel[[which(unique(fam)==f)]])
            ww[fam==f,,,,] <- res$ww
            w[fam==f,,] <- res$w
            ll <- ll+res$ll
        }
        res <- list("ww"=ww,"w"=w,"ll"=ll)
    }
    else
    {
        for(i in 1:length(id))
        {
            if(status[i]==2)
            {
                w[i,1,] <- fyc[i,]*c(probs$p,0)*probs$p.aff
                w[i,2,] <- 0
                ll <- ll+log(sum(w[i,,]))
            }
            if(status[i]==1)
            {
                w[i,1,] <- 0
                w[i,2,] <- c((1-probs$p0)*probs$p,probs$p0)*(1-probs$p.aff)
                ll <- ll+log(sum(w[i,,]))
            }
            if(status[i]==0)
            {
                w[i,1,] <- c(probs$p,0)*probs$p.aff
                w[i,2,] <- c((1-probs$p0)*probs$p,probs$p0)*(1-probs$p.aff)
                ll <- ll+log(sum(w[i,,]))
            }
            w[i,,] <- w[i,,]/sum(w[i,,])
        }
        res <- list("w"=w,"ll"=ll)
   }
   cat("log-likelihood : ",ll,"\n")
   res
}

