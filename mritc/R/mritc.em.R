emnormmix <- function(y, prop, mu, sigma, err, maxit, verbose){
    checkErrors(prop, mu, sigma, err)

    if (length(err) < 3) err <- rep(err, length.out = 3)

    k <- length(prop)
    ## tabulate distinct y values
    yvals <- sort(unique(y))
    n.yvals <- length(yvals)
    counts <- as.vector(table(y))

    ## enlarge y to dy.
    dy <- rep(yvals,each=k)
    no <- length(y)

    it <- 0

    repeat {
        it <- it + 1
        muold <- mu
        sigmaold <- sigma
        propold <- prop

        ## compute new component probabilities
        dprop <- rep(propold, n.yvals)
        dmu <- rep(muold, n.yvals)
        dsigma <- rep(sigmaold, n.yvals)
        dn <- dnorm(dy, mean = dmu, sd = dsigma)
        interprop <- matrix(dprop * dn, ncol = k, byrow = TRUE)
        indices <- sweep(interprop, 1, rowSums(interprop), `/`)

        ## compute new weights
        sumIndices <- as.vector(counts %*% indices)
        prop <- sumIndices/no

        ## compute new means
        mu <- as.vector(((yvals * counts) %*% indices)/sumIndices)

        ## compute new variances
        sigma <- sqrt(counts %*% (outer(yvals, mu, `-`) ^ 2 * indices) / sumIndices)

        flag <- checkStopVerbose(muold, mu, sigmaold, sigma, err, it, maxit, verbose, propold, prop)
        if(flag==1) break

      
    }
    list(prop = prop, mu = mu, sigma = sigma)
}


#############################################################
mritc.em <- function(y, prop, mu, sigma, err=1e-4, maxit=200, verbose) {

    checkErrors(prop, mu, sigma, err)

    k <- length(prop)
    nvert <- length(y)

    #get the EM estimation
    result <- emnormmix(y, prop, mu, sigma, err, maxit, verbose)
    prop <- result$prop
    mu <- result$mu
    sigma <- result$sigma


    ## tabulate distinct y values
    yvals <- sort(unique(y))
    n.yvals <- length(yvals)
    counts <- as.vector(table(y))

    ## k : number of normal distributions in the mixture.
    k <- length(prop)

    ## enlarge y to dy.
    dy <- rep(yvals,each=k)
    no <- length(y)

    ## compute component probabilities
    dprop <- rep(prop, n.yvals)
    dmu <- rep(mu, n.yvals)
    dsigma <- rep(sigma, n.yvals)
    dn <- dnorm(dy, mean = dmu, sd = dsigma)
    interprop <- matrix(dprop * dn, ncol = k, byrow = TRUE)
    indices <- sweep(interprop, 1, rowSums(interprop), `/`)

    pre <- indices

    list(prob=pre[match(y, yvals),], mu=mu, sigma=sigma)
}
