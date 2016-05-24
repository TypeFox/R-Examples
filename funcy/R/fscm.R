#
#  Copyright (C) 2011-2015 Christina Yassouridis
#  
#

fscm <- function(data, k, reg, regTime, funcyCtrlMbc, fpcCtrl,
                 location, scale, knn, useCode, verbose){

    ##control parameters
    dimBase <- funcyCtrlMbc@dimBase
    baseType <- funcyCtrlMbc@baseType
    eps <- funcyCtrlMbc@eps
    maxit <- funcyCtrlMbc@maxit
    seed <- funcyCtrlMbc@seed
    init <- funcyCtrlMbc@init
    nrep <- funcyCtrlMbc@nrep
 
    nt  <-  dim(data)[1]
    nc  <-  dim(data)[2]

    if(is.null(location)){
        gridNr <- squareGrid(x=nc)
        location <- expand.grid(1:gridNr[1], 1:gridNr[2])
    }
    
    ##scale dataset
    if(scale)
        data <- apply(as.matrix(data), 2, function(x) x-mean(x))
    data <- t(data)
   
    ##calculate spatial covariance matrix
    B <- stationary.cov(location,  location,  Covariance="Matern")
    psi <- inchol(B)
    psi <- t(psi)
   
    ##calculate basis functions
    if(is.null(regTime))
            grid <- 1:nt
        else
            grid <- regTime

    ##phi matrix
    if(baseType=="eigenbasis"){
        mkdat <- formatFuncy(data=t(data), regTime=regTime,
                             format="Format3", reg=TRUE)
        res <- fpc(Yin=mkdat$Yin, Tin=mkdat$Tin, isobs=mkdat$isobs, dimBase=dimBase,
                   fpcCtrl=fpcCtrl, reg=TRUE)
        plotParams <- res$plotParams
        phi <- res$base
    }else{
        res <- makeBasis(basis=baseType, time=grid, nbasis=dimBase)
        bObj <- res$bObj
        phi <- res$phi
    }
    
    phi <- svd(phi)$u
    phi <- t(phi)

    ##obtain neighborhood matrix: input 'nb'
    dist <- as.matrix(dist(location))
    dist.ord <- apply(dist, 2, order)
    dist.order <- t(dist.ord)
    nb <- dist.ord[, 2:(knn+1)]-1
    nb <- t(nb)

    ##initial clustering
    cl <- initClust(data, k, init=init, seed=seed, nrep=nrep)$cluster-1
    z.init <- t(cl)

    if(useCode=="R"){
        z.init <- as.integer(z.init+1)
        nb <- t(nb)
        res=Rfscm((data), nb, t(phi), t(psi), z.init, maxit)
        k <- res$k
        alpha <- rowMeans(sapply(1:nc, function(j)
            ginv(phi%*%t(phi))%*%(phi%*%(data[j,]-t(phi)%*%res$beta[,res$z[j]]-rep(t(psi[j,])%*%res$gamma, nt)))))
        overall <- as.numeric(t(phi)%*%alpha)
        temp <- t(phi)%*%res$beta
        centers <- overall + temp
        temp <- cbind(overall, temp)
        spatial <- as.numeric(t(psi)%*%res$gamma)
        trends <- list(temp=temp,  spatial=spatial)
        probs <- res$prob_z
        cls <- res$z
        params <- list(beta=res$beta, gamma=res$gamma,
                       sigma_e=res$sigma_e, sigma_s=res$sigma_s,
                       theta=res$theta, trends=trends)
        AIC <- BIC <- NA
    
    }else if(useCode=="C"){
        temps <- tempfile()
        dir.create(temps)
        
        ##create a temporary directory
        dataFile <- paste0(temps, "/data.txt")
        write.table((data),dataFile,row.names=F, col.names=F)
        phiFile <- paste0(temps, "/phi.txt")
        write.table(phi,  file=phiFile, row.names=F, col.names=F)
        psiFile <- paste0(temps, "/psi.txt")
        write.table(psi, psiFile, row.names=F, col.names=F)
        nbFile <- paste0(temps, "/nb.txt")
        write.table(nb,  nbFile, row.names=F, col.names=F)
        zFile <- paste0(temps, "/z.init", k, ".txt")
        write.table(z.init, zFile, row.names=F, col.names=F)


        olddir <- getwd()
        setwd(temps)
        on.exit(setwd(olddir))
       
        res <- .C("do_fscm", as.character(dataFile), as.character(nbFile), 
            as.character(zFile), as.character(phiFile), as.character(psiFile),
            as.integer(k), as.double(eps), 
            1000L, # max_iter_gamma
            100L,  # max_iter_z
            as.integer(maxit), 
            100L,  # max_iter_hmrf
            as.integer(verbose),
            ok=0L)  #error    

        # Error occured during the estimation.
        if(!res$ok) {
            stop("numerical problems in fscm, restart the computation with a different seed")
        }

        AIC <- BIC <- NA
        trends <- params <- list()
        tp <- read.csv2("temporal.csv", sep="\t",  header=FALSE)
        trends$temp <- do.call(rbind, lapply(tp, function(x) as.numeric(gsub("[,;]", "", x))))
        trends$spatial <- read.csv("spatial.csv", header=FALSE,sep=",")
        params$alpha <- read.table("alpha.csv", header=FALSE, sep=",")[,-1]
        params$beta <- read.table("beta.csv", header=FALSE, sep=",")
        params$gamma <- read.table("gamma.csv", header=FALSE,
                                   sep=",")[,-1]
        params$sigma <- trends$spatial[,1]
        params$sigmaErrors <- t(read.csv("parameter.csv", header=FALSE,
                                         sep=",")[,-1])
        colnames(params$sigmaErrors) <- c("sigma_e","simga_s")
        params$theta <- as.numeric(read.table("theta.csv",
                                              header=FALSE))

        setwd(olddir)
  
        cls <- as.numeric(trends$spatial[2,])+1
        centers <- trends$temp[,-1]+trends$temp[,1]
    }
    
    dev <- rep(0,nc)
    for(i in 1:k) {
        ind <- which(cls==i)
        ctrs <- trends$temp[,i]
        ctrs <- centers[,i]
        dat <- data[ind,]
        diff <- t(t(dat)-ctrs)
        if(length(ind)==1)
            diff <- t(diff)
        dev[ind] <- apply(diff,1,function(x) sum(x^2))
    }
    params$sigma <- (dev-min(dev))/(max(dev)-min(dev))

    isNA <- which(apply(centers, 2, function(x) sum(is.na(x)))==nt)
    if (length(isNA) != 0) {
        k <- k - length(isNA)
        warning("Class number was reduced to ",k,". Try with another seed or useCode=R.")
        centers <- centers[, -isNA, drop = FALSE]
        params$beta <- params$beta[-isNA, , drop = FALSE]
        trends$temp <- trends$temp[, -(isNA + 1), drop = FALSE]
    }
    
    res <- list(methodName="FSCM",prms=params,k=k,dimBase =
                    dimBase,cluster=cls,centers=centers, data =
                        data, AIC=AIC, BIC=BIC, location=location,
                trends=trends, grid=grid)
    
    return(res) 
}


