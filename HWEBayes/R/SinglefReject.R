SinglefReject <-
function(nsim,bvec,lambdamu,lambdasd,nvec){
    k <- length(bvec)
    if (length(nvec) != k*(k + 1)/2) {
        stop("length mismatch between bvec and nvec")
    }
    MLEres <- HWEmodelsMLE(nvec)
    maxLL <- MLEres$fmaxloglik

    ##psamp <- matrix(0,nrow=nsim,ncol=k)
    ##fsamp <- NULL
    ##psamp <- fsamp <- vector("list", nsim)

    psamp <- fsamp <- vector("list", nsim)

    niter <- naccept <- count <- PrnH1 <- varterm1 <- 0
    Lnorm <- lfactorial(sum(nvec)) - sum(lfactorial(nvec))

    while (naccept < nsim){

        M <- nsim - naccept
        count <- count + M
        samples <- SinglefPrior(nsim=M,alpha=bvec,lambdamu=lambdamu,
                                lambdasd=lambdasd)

        LL <- if(M==1) {
            MultLogLikP(samples$p, samples$f, nvec)
        }
        else {
            apply(cbind(samples$f, samples$p), 1,
                  function(x) MultLogLikP(x[-1], x[1], nvec))
        }
        
        likterm <- LL + Lnorm
        expterm <- exp(likterm)
        PrnH1 <- PrnH1 + sum(expterm)
        varterm1 <- varterm1 + sum(expterm^2)
        if (any(LL > maxLL)) cat("Maximization is messed up\n")
        accept <- log(runif(M)) < LL - maxLL
        if (any(accept)) {
            niter <- niter + 1
            naccept <- naccept + sum(accept)
            psamp[[niter]] <- samples$p[accept,,drop=FALSE]
            fsamp[[niter]] <- samples$f[accept]
        }
###if (floor(naccept/100)==ceiling(naccept/100))cat("Number of accepted points = ",naccept,"\n")
    }
    psamp <- do.call(rbind,psamp)
    fsamp <- unlist(fsamp)

    PrnH1 <- PrnH1/count
    varest <- (varterm1/count - PrnH1^2)/count
    cat("nsim norm constant (se) 95% interval: \n")
    cat(nsim,PrnH1,"(",sqrt(varest),")",PrnH1-1.96*sqrt(varest),PrnH1+1.96*sqrt(varest),"\n")
    accrate <- nsim/count
    list(psamp=psamp,fsamp=fsamp,accrate=accrate,PrnH1=PrnH1,varest=varest)
}

