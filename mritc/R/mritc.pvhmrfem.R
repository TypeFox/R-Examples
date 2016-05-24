inte <- function(yunique, n.yunique, ymatch, mu, sigma, n){
    n.yunique <- length(yunique)
    alpha <- seq(0, 1, length=n+1)
    sigma2 <- sigma^2
    mu <- alpha * mu[1] + (1 - alpha) * mu[2]
    sigma <- sqrt(alpha^2 * sigma2[1] + (1 - alpha)^2 * sigma2[2])
    dn <- getDen(yunique, n.yunique, 1:n.yunique, mu, sigma)
    den <- apply(dn, 1, function(x) (sum(2 * x) - x[1] - x[n+1]) / (2 * n) )
    den <- den[ymatch]
    den

}

mritc.pvhmrfem <- function(y, neighbors, blocks,
                      spatialMat=matrix(c(2, 1, -1, -1, -1, 1, 2, 1, -1, -1,
                        -1, 1, 2, 1, -1, -1, -1, 1, 2, 1, -1, -1, -1, 1, 2),
                        ncol=5), beta=0.6,
                      mu, sigma, err=1e-4, maxit=20, verbose){
    #Error checking

    if (length(err) < 2) err <- rep(err, length.out = 2)

    k <- length(mu) + 2
    nneigh <- ncol(neighbors)
    nblocks <- length(blocks)
    neighbors <- structure(as.integer(neighbors), dim = dim(neighbors))

    yunique <- sort(unique(y))
    n.yunique <- length(yunique)
    nvert <- length(y)
    ymatch <- match(y, yunique)


    muold <- rep(0,k-2)
    sigmaold <- rep(0,k-2)
    niter <- 10

    check <- getCheck(nneigh, k, beta, spatialMat)
    indices <- initialIndices(y, nvert, mu, sigma, k, sub=FALSE, type="partial")

    it <- 0
    repeat{
        it <- it + 1
        muold <- mu
        sigmaold <- sigma
        den.part <- getDen(yunique, n.yunique, ymatch, mu, sigma)
        den <- matrix(0, nrow=nvert, ncol=k)
        den[,1] <- den.part[,1]
        den[,3] <- den.part[,2]
        den[,5] <- den.part[,3]
        den[,2] <- inte(yunique, n.yunique, ymatch, mu[1:2], sigma[1:2], 100)
        den[,4] <- inte(yunique, n.yunique, ymatch, mu[2:3], sigma[2:3], 100)

        for(j in 1:niter){
            indices <- updateIndicesHMRFEM(neighbors, nneigh, blocks, nblocks, as.integer(k), indices, check, den)
        }
        prob  <- updateProbs(neighbors, nneigh, indices, den, k, check)
        mu <- updateMus(prob[,c(1,3,5)], y)
        sigma <- updateSds(prob[,c(1,3,5)], y, mu, nvert, k-2)
        flag <- checkStopVerbose(muold, mu, sigmaold, sigma, err, it, maxit, verbose)
        if(flag==1) break
    }

    list(prob=prob, mu=mu, sigma=sigma)
}

