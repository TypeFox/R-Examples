#get possible configurations of a graph
getConfs <- function(n, nc){
    if(n < 2) stop("There are at least 2 vertices!")
    if(nc < 2) stop("There are at least two possible choices for each vertex!")

    subcc <- function(num){
        lapply(1:nc, function(x) list(c(num, x)))
    }

    cc <- function(n){
        if(n > 2){
            part <- matrix(unlist(cc(n-1)), nrow=n-1)
            lapply(1:ncol(part), function(x) subcc(part[,x]))
        }
        else{
            lapply(1:nc, function(x) subcc(x))
        }

    }

    matrix(unlist(cc(n)), nrow=n)
}
#get unique configurations 
getConfsUni <- function(nneigh, k){
    com <- getConfs(k, nneigh+1)
    com <- com - 1
    or <- order(apply(com, 2, function(x)
                      sum(sapply(1:k, function(i) x[i]*(nneigh+1)^(i-1)))))
    com[,or]
}

#compute the check table of probabilities for different neighbor configurations used in updating the indices based on a Potts model
getCheck <- function(nneigh, k, beta, spatialMat){
    if(k==3){
        if(nneigh==6) load(system.file("tables/tabn6k3.rda", package="mritc"))
        else
            if(nneigh==18) load(system.file("tables/tabn18k3.rda", package="mritc"))
            else
                if(nneigh==26) load(system.file("tables/tabn26k3.rda", package="mritc"))
                else com <- getConfsUni(nneigh, k)     
    }
    else{
        com <- getConfsUni(nneigh, k) 
    }
    exp(beta *   t(apply(com, 2, function(x) spatialMat%*%x)))
}

#check errors of input of tissue classi 
checkErrors <- function(prop=NULL, mu, sigma, err=NULL){
    if(is.null(prop)){
        k <- length(mu)
        if (k != length(sigma))
            stop("The dimensions of 'mu' and 'sigma' do NOT match.")
        if (!all(sigma > 0))
            stop("All 'sigma's have to be positive")
    }
    else{
        if (!(all(prop > 0) && all(prop < 1)))
            stop("'prop' has to be between 0 and 1")
        if (!isTRUE(all.equal(sum(prop), 1)))
            stop("Sum of 'prop' has to be 1")
        if (!is.null(mu) && !is.null(sigma)){
            k <- length(prop)
            if (!(k == length(mu) && k == length(sigma)))
                stop("The dimensions of 'prop', 'mu' and 'sigma' do NOT match.")
        }
    }
    if (!all(sigma > 0))
        stop("All 'sigma's have to be positive")
    if(! is.null(err))
        if (any(err < 0)) stop("All 'err's have to be positive.")

}

#Compute the density of all vertices.
getDen <- function(yunique, n.yunique, ymatch, mu, sigma){
    k <- length(mu)
    dy <- rep(yunique, each=k)
    dmu <- rep(mu, n.yunique)
    dsigma <- rep(sigma, n.yunique)
    den <-  dnorm(dy, mean = dmu, sd = dsigma)
    den <- matrix(den, ncol=k, byrow=T)
    den <- den[ymatch,]
    den
}

#compute the relative errors
relerr <- function(x, y){
    max(abs(x - y)) / (1 + max(abs(x), abs(y)))
}

#Initialize indices
initialIndices <- function(y, nvert, mu, sigma, k, sub, subvox=NULL, type="pure"){
    if(sub == FALSE){
        if(type == "pure")
            indices <- max.col(matrix(dnorm(rep(y, k), mean=rep(mu, each=nvert),
                                            sd=rep(sigma, each=nvert)), ncol=k))
        else{
            indices <- matrix(dnorm(rep(y, k-2),mean=rep(mu, each=nvert),
                          sd=rep(sigma, each=nvert)), ncol=k-2)
            indices <- round((indices[,1] + 3*indices[,2] + 5*indices[,3]) /rowSums(indices))
        }

    }
    else{
        nobs <- nvert / 8
        index <- max.col(matrix(dnorm(rep(y/8,k),
                                      mean=rep(mu, each=nobs),
                                      sd=rep(sigma, each=nobs)), ncol=k))
        indices <- rep(0, nvert)
        indices[subvox] <- rep(index, 8)
    }
    indices <- do.call(cbind, lapply(1:k, function(x) indices==x))
    rbind(indices, rep(0L,k))
 
}
 
#check whether to stop the interation or not and whether ouput the curren state
checkStopVerbose <- function(muold, mu, sigmaold, sigma, err, it, maxit,
                             verbose, propold=NULL, prop=NULL){
    flag <- 0
    remu <- relerr(muold, mu)
    resigma <- relerr(sigmaold, sigma)
    if(is.null(prop)){
        if (remu < err[1] && resigma < err[2] || it > maxit)
            flag <- 1
        if (verbose)
            cat(paste("Iteration ", it,
                      ": relative error of mu = ", signif(remu, 1),
                      "; sigma = ", signif(resigma, 1),
                      "\n", sep=""))
    }
    else{
        reprop <- relerr(propold, prop)
        if (remu < err[1] && resigma < err[2] && reprop < err[3] || it > maxit)
            flag <- 1
        if (verbose && it %% 10 == 0)
            cat(paste("Iteration ", it,
                      ": relative error of mu = ", signif(remu, 1),
                      "; sigma = ", signif(resigma, 1),
                      "; prop = ", signif(reprop, 1),
                      "\n", sep=""))
        
    }
    flag
}


    
