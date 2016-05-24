mritc.hmrfem <- function(y, neighbors, blocks, spatialMat=diag(1,3), beta=0.5, mu, sigma, err=1e-4, maxit=20, verbose){

    checkErrors(mu=mu, sigma=sigma, err=err)

    if (length(err) < 2) err <- rep(err, length.out = 2)

    k <- length(mu)
    nvert <- length(y)


    nneigh <- ncol(neighbors)
    nblocks <- length(blocks)
    neighbors <- structure(as.integer(neighbors), dim = dim(neighbors))

    yunique <- sort(unique(y))
    n.yunique <- length(yunique)
    nvert <- length(y)
    ymatch <- match(y, yunique)

    muold <- rep(0,k)
    sigmaold <- rep(0,k)
    niter <- 10

    check <- getCheck(nneigh, k, beta, spatialMat)
    indices <- initialIndices(y, nvert, mu, sigma, k, sub=FALSE)

    it <- 0

    repeat{
        it <- it + 1
        muold <- mu
        sigmaold <- sigma
        for(j in 1:niter){
            den <- getDen(yunique, n.yunique, ymatch, mu, sigma)
            indices <- updateIndicesHMRFEM(neighbors, nneigh, blocks, nblocks, k, indices, check, den)
        }
        prob  <- updateProbs(neighbors, nneigh, indices, den, k, check)
        mu <- updateMus(prob, y)
        sigma <- updateSds(prob, y, mu, nvert, k)

        flag <- checkStopVerbose(muold, mu, sigmaold, sigma, err, it, maxit, verbose)
        if(flag==1) break
    }

    list(prob=prob, mu=mu, sigma=sigma)
}

