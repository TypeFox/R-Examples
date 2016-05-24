initalBias <- function(y, indices, mu, nvert, blocks, neighbors, nneighbors, weights){
    bias <- y - (indices %*% mu)[-(nvert+1)]
    bias <- c(bias,0)
    for(i in 1:length(blocks)){
        points <- blocks[[i]]
        mat <- matrix(bias[neighbors[points,]], nrow=length(points))
        mat <- mat %*% weights
        bias[points] <- (mat +  bias[points]) /
          (nneighbors[points] + 1)
    }

    bias

}

MCMCSubclass <- function(indices, niter, subvox){
    ind <- indices[t(subvox),]
    nmask <- nrow(subvox)
    p <- matrix(0, nrow=nmask,ncol=ncol(indices))
    for(i in 1:8){
        a <- seq(1,(nmask*8)-7,by=8) + i - 1
        p <- p + ind[a,]
    }
    p.vox <- p/8/niter

    ind <- indices[t(subvox),]
    mode <- matrix(max.col(ind),nrow=8)
    class.modemean=cbind(colSums(mode==1),colSums(mode==2), colSums(mode==3))/8

    list(fussy.vox = p.vox, class.modemean=class.modemean)
}

updateY <- function(yobs, indices, subvox, mu, sigma)
    .Call("updateY", yobs, subvox, indices, mu, sigma)

updateYBias <- function(yobs, bias, indices, subvox, mu, sigma)
    .Call("updateYBias", yobs, bias, subvox, indices, mu, sigma)

getDenSub <- function(y, mu, sigma)
    .Call("getDenSub", y, mu, sigma)

getDenSubBias <- function(y, bias, mu, sigma)
    .Call("getDenSubBias", y, bias, mu, sigma)

updateIndices <- function(neighbors, nneigh, blocks, nblocks, k, indices, check,den){
    .Call("updateIndices", blocks, neighbors, nneigh, k, indices, check, den)
}

updateMusB <- function(Nj, S1j, sigma, eta, xi, k){
    oldtau <- 1/(sigma^2)
    newpre <- eta + oldtau * Nj
    newmean <- (eta*xi + oldtau * S1j) / newpre
    mu <- rnorm(k, mean=newmean, sd=1/sqrt(newpre))
    mu
}

updateSdsB <- function(Nj, ybar, S2j, mu, lambda, phi, k){
    shape <- Nj/2 + lambda
    mid <- S2j + Nj * (ybar - mu)^2
    scale <- 2 * phi / (2 + phi*mid)
    newtau <- rgamma(k, shape=shape, scale=scale)
    sigma <- 1/sqrt(newtau)
    sigma
}

yIndicesSummaries <- function(y, indices, k)
    structure(.Call("yIndicesSummaries", y, indices, k), names = c("Nj", "ybar", "S2j"))


yIndicesSummariesBias <- function(y, bias, indices, k)
    structure(.Call("yIndicesSummariesBias", y, bias, indices, k), names = c("Nj", "ybar", "S2j"))



updateSmoothWeights <- function(gamma.a, gamma.b, nedges,
                                neighbors, nneigh, nvertex,
                                bias, weights){
    .Call("updateSmoothWeights", gamma.a, gamma.b, nedges,
          neighbors, nneigh, nvertex, bias, weights)

}

updateBiasInde <- function(y, bias, neighbors, nneigh, weineighbors, indices,
                       mu, sigma, gamma, weights){
    .Call("updateBiasInde", y, bias, neighbors, nneigh, weineighbors, indices,
          mu, sigma, gamma, weights)
}

updateBias <- function(y, bias, blocks, neighbors, nneigh, weineighbors, indices,
                       mu, sigma, gamma, weights){
    .Call("updateBias", y, bias, blocks, neighbors, nneigh, weineighbors, indices,
          mu, sigma, gamma, weights)
}

updateDistanceUnit <- function(ubsigma.omega, nneigh, nvert, gamma, sigma.omega, sdomega, bias, neighbors, weights, weineighbors, neiDiscrete){

    sigma.omega.new <- qnorm(runif(1, pnorm(0, sigma.omega, sdomega),
                                    pnorm(ubsigma.omega, sigma.omega, sdomega)),
                              sigma.omega, sdomega)
    weights.new <- getWeightsMRI(nnei=nneigh, sigma=sigma.omega.new)
    weineighbors.new <-  rowSums(neiDiscrete %*% weights.new)

    qq <- .Call("updateDistanceUnit", bias, neighbors, nvert, nneigh, weights, weights.new, weineighbors, weineighbors.new, neiDiscrete)
    ratio <- gamma/2*(qq[1] - qq[2]) +
             log(pnorm(ubsigma.omega, sigma.omega, sdomega)-pnorm(0, sigma.omega, sdomega)) -
             log(pnorm(ubsigma.omega, sigma.omega.new, sdomega)-pnorm(0, sigma.omega.new, sdomega))
    if(! is.na(ratio)){
        if(ratio > log(runif(1))){
            sigma.omega <- sigma.omega.new
            weights <- weights.new
            weineighbors <- weineighbors.new
        }
    }
    list(sigma.omega=sigma.omega, weights=weights, weineighbors=weineighbors)


}

mritc.bayes.nobias <- function(y, neighbors, blocks, sub, subvox,
                        spatialMat, beta, mu, sigma, niter, verbose){
    #r <- range(yobs)
    r <- 256
    xi <- r/2
    m1 <- 0.3
    eta <- m1/r^2
    lambda <- m1
    phi <- 1/r^2

    checkErrors(mu=mu, sigma=sigma, err=NULL)

    k <- length(mu)

    nneigh <- ncol(neighbors)
    nblocks <- length(blocks)
    neighbors <- structure(as.integer(neighbors), dim = dim(neighbors))

    y <- as.double(y)

    check <- matrix(as.double(getCheck(nneigh, k, beta, spatialMat)), ncol=k)

    if(sub==F){
        yunique <- sort(unique(y))
        n.yunique <- length(yunique)
        nvert <- length(y)
        ymatch <- match(y, yunique)
        indices <- initialIndices(y, nvert, mu, sigma, k, sub=FALSE)
    }

    else{
        yobs <- y
        nobs <- length(yobs)
        nvert <- nobs*8
        mu <- mu/8
        sigma <- sqrt(sigma^2)/8
        y <- rep(0, nvert)
        y[subvox] <- rep(yobs/8, 8)
        xi <- xi/8
        eta <- eta*64
        lambda <- lambda
        phi <- phi*64
        indices <- initialIndices(yobs, nvert, mu, sigma, k, sub=TRUE, subvox)
    }


    indicesMCMC <- matrix(0L, nrow=nrow(indices), ncol=ncol(indices))
    musave <- matrix(0, ncol=k, nrow=niter)
    sigmasave <- matrix(0, ncol=k, nrow=niter)

    for (i in 1:niter){
        if(sub==F)
            den <- getDen(yunique, n.yunique, ymatch, mu, sigma)
        else den <- getDenSub(y, mu, sigma)

        indices <- updateIndices(neighbors, nneigh, blocks, nblocks, k, indices,
                           check, den)

        indicesMCMC <- .Call("updateCounts", indicesMCMC, indices)

        if(sub==T){
            y <- updateY(yobs, indices, subvox, mu, sigma)
        }

        yIndicesSums <- yIndicesSummaries(y, indices, k)
        Nj <- yIndicesSums$Nj
        ybar <- yIndicesSums$ybar
        S2j <- yIndicesSums$S2j
        S1j <- Nj * ybar
        mu <- updateMusB(Nj, S1j, sigma, eta, xi, k)
        sigma <- updateSdsB(Nj, ybar, S2j, mu, lambda, phi, k)

        musave[i,] <- mu
        sigmasave[i,] <- sigma

                
        if(verbose){
            percent <- i / (niter) * 10
            if(percent %% 1 == 0)
              cat(paste(percent*10, "%", " of iterations has finished", "\n", sep=""))
        }

    }
    if(sub==FALSE)
        prob <- indicesMCMC[-(nvert+1),]/niter
    else
        prob <- class <- MCMCSubclass(indicesMCMC[-(nvert+1),], niter, subvox)$fussy.vox
    list(prob=prob, mu=colMeans(musave), sigma=colMeans(sigmasave))
}


mritc.bayes.bias <- function(y, neighbors, blocks, subvox,
                             neighbors.bias, blocks.bias,
                             weineighbors.bias, weights.bias,
                             spatialMat=matrix(c(2,0,-1,0,2,0,-1,0,2), nrow=3),
                             beta=0.3, mu, sigma, niter=1000, verbose=TRUE){
    r <- 256
    xi <- r/2
    m1 <- 0.3
    eta <- m1/r^2
    lambda <- m1
    phi <- 1/r^2
    gamma.a <- gamma.b <- 1
    nInde <- 200
    ubsigma.omega <- 10
    sigma.omega <- 1
    sdomega <- 0.1
    
    checkErrors(mu=mu, sigma=sigma, err=NULL)
    
    k <- length(mu)
    nneigh <- ncol(neighbors)
    nneigh.bias <- ncol(neighbors.bias)
    nblocks <- length(blocks)
    neighbors <- structure(as.integer(neighbors), dim = dim(neighbors))
    neighbors.bias <- structure(as.integer(neighbors.bias), dim = dim(neighbors.bias))
    weineighbors.bias <- structure(as.double(weineighbors.bias), dim = dim(weineighbors.bias))
    nedges.bias <- as.integer(sum(neighbors.bias != nrow(neighbors.bias)+1) / 2)
    
    y <- as.double(y)
    check <- matrix(as.double(getCheck(nneigh, k, beta, spatialMat)), ncol=k)
    yobs <- y
    nobs <- length(yobs)
    nvert <- as.integer(nobs*8)
    mu <- mu/8
    sigma <- sqrt(sigma^2)/8
    y <- rep(0, nvert)
    y[subvox] <- rep(yobs/8, 8)
    xi <- xi/8
    eta <- eta*64
    lambda <- lambda
    phi <- phi*64
    neiDiscrete <- (neighbors.bias < nvert+1)
    neiDiscrete <- structure(as.integer(neiDiscrete), dim=dim(neiDiscrete))
    indices <- initialIndices(yobs, nvert, mu, sigma, k, sub=TRUE, subvox)
    bias <- initalBias(y, indices, mu, nvert, blocks.bias, neighbors.bias, weineighbors.bias, weights.bias)
    
    
    indicesMCMC <- matrix(0L, nrow=nrow(indices), ncol=ncol(indices))
    musave <- matrix(0, ncol=k, nrow=niter)
    sigmasave <- matrix(0, ncol=k, nrow=niter)
    
    for (i in 1:(niter+nInde)){
        den <- getDenSubBias(y, bias, mu, sigma)

        indices <- updateIndices(neighbors, nneigh, blocks, nblocks, k, indices,
                                 check, den)
    
        if(i > nInde){
            indicesMCMC <- .Call("updateCounts", indicesMCMC, indices)
        }
    
        y <- updateYBias(yobs, bias, indices, subvox, mu, sigma)
    
        yIndicesSums <- yIndicesSummariesBias(y, bias, indices, k)
        Nj <- yIndicesSums$Nj
        ybar <- yIndicesSums$ybar
        S2j <- yIndicesSums$S2j
        S1j <- Nj * ybar
        mu <- updateMusB(Nj, S1j, sigma, eta, xi, k)
        sigma <- updateSdsB(Nj, ybar, S2j, mu, lambda, phi, k)
        
        if(i > nInde){
            musave[(i-nInde),] <- mu
            sigmasave[(i-nInde),] <- sigma
        }
        
        gamma <- updateSmoothWeights(gamma.a, gamma.b, nedges.bias,
                                     neighbors.bias,
                                     nneigh.bias, nvert,
                                     bias, weights.bias)
        #if(i > nInde){
        #    gammasave[(i-nInde)] <- gamma
        #}
        
        if(i <= nInde){
            bias <- updateBiasInde(y, bias, neighbors.bias, nneigh.bias,
                                   weineighbors.bias, indices,
                                   mu, sigma, gamma, weights.bias)
        }
        else{
            bias <- updateBias(y, bias, blocks.bias, neighbors.bias, nneigh.bias,
                               weineighbors.bias, indices, mu, sigma, gamma,
                               weights.bias)
        }
        
        if(i > nInde){
            sigma.omega.old <- sigma.omega
            omega <- updateDistanceUnit(ubsigma.omega, nneigh.bias, nvert,
                                        gamma, sigma.omega,
                                        sdomega, bias, neighbors.bias,
                                        weights.bias, weineighbors.bias,
                                        neiDiscrete)
            sigma.omega <- omega$sigma.omega
            weights.bias <- omega$weights
            weineighbors.bias <- omega$weineighbors
            
            
            #sigmaomegasave[(i-nInde)] <- sigma.omega
        }
        
        if(verbose){
            percent <- i / (niter+nInde) * 10
            if(percent %% 1 == 0)
              cat(paste(percent*10, "%", " of iterations has finished", "\n", sep=""))
        }
    }
    
    prob <- class <- MCMCSubclass(indicesMCMC[-(nvert+1),], niter, subvox)$fussy.vox
    list(prob=prob, mu=colMeans(musave), sigma=colMeans(sigmasave))
}

mritc.bayes <- function(y, neighbors, blocks, sub, subvox,
                        subbias=FALSE, neighbors.bias=NULL, blocks.bias=NULL,
                        weineighbors.bias=NULL, weights.bias=NULL,
                        spatialMat=(if(sub) matrix(c(2,0,-1,0,2,0,-1,0,2), nrow=3)
                                    else diag(1,3)),
                        beta=ifelse(sub, 0.3, 0.7), mu, sigma, niter=100, verbose=TRUE){
    if(subbias==FALSE){
        mritc.bayes.nobias(y, neighbors, blocks, sub, subvox,
                           spatialMat, beta, mu, sigma, niter, verbose)
    }
    else{
        if(sub==FALSE){
            stop("PV and intensity non-uniformity have to be addressed simultaneously!")
        }
        else{
            mritc.bayes.bias(y, neighbors, blocks, subvox,
                             neighbors.bias, blocks.bias,
                             weineighbors.bias, weights.bias,
                             spatialMat, beta, mu, sigma, niter, verbose)
        }
    }
    
}
