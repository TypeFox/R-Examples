### Original is  /usr/local/app/R/R_local/src/covRobust/R/cov.nnve.R
###              9129   Nov 18 2003  cov.nnve.R
### --> keep as covNNC1 (unexported) in ./covNNC-orig.R

## Much improved by Martin Maechler :

##- speed:    Mtovec, vectoM sped up; later eliminated
##            to do the *scaling* once for X, then for the n2 rows
##- accuracy: Ddk(*, lambda1,*) & Ddk(*, lambda2, *) get very close to 0
##            But we can  use q.Ddk() := Ddk(*, ..1) / Ddk(*, ..2)  with
##            much more accurate formula.
## - Fix to make it work for d = ncol(X) = 1

## Simply renamed   (for readability / understandability / ...):
## s/datamat/X/
## s/SplusNN/NNdist/

covNNC <- function(X, k = min(12, n-1), pnoise = 0.05,
                   emconv = 0.001, bound = 1.5, extension = TRUE, devsm = 0.01)
{
    ## Function to perform Nearest Neighbor Variance Estimation
    ## Wang and Raftery(2002), "Nearest neighbor variance estimation (NNVE):
    ##     Robust covariance estimation via nearest neighbor cleaning
    ##                 (with discussion)",
    ##     Journal of the American Statistical Association 97:994-1019
    ## Available as Technical Report 368 (2000) from
    ##    http://www.stat.washington.edu/www/research/reports
    ##
    ## The data need to be a matrix of points in the process.
    ## Each column contains one variable.
    ##
    ##
    ##  Returns a list with the following components:
    ##      cov            - robust covariance estimate
    ##      mu             - mean
    ##      postprob       - posterior probability
    ##      classification - classification (0=noise otherwise 1)
    ##                      (obtained by rounding postprob)
    ##      innc - list of initial nearest-neighbor results
    ##            (components are the same as above)
    ##
    ##---------------------------------------------------------------------------

    ##
    ##  The Re-scale NNC function
    ##
    nclean.sub <- function(X, k, convergence = emconv,
                           s.X, # := the centered & scaled X
                           k.dists = NNdist(s.X, k = k))
    {
        ##
        ## Function to perform the Nearest Neighbour cleaning of
        ##--------------------------------------------------------------

        dDk <- function(x, lambda, k, d, alpha.d) {
            ## find the density of D_k
            exp( - lambda * alpha.d * x^d + log(2) +
                k * log(lambda * alpha.d) + log(x) * (d * k - 1) - log(gamma(k)))
        }

        ## MM:
        q.dDk <- function(x, lam1, lam2, k, d, alpha.d) {
            ## := dDk(x, lam1, k,d,alpha.d) /
            ##    dDK(x, lam2, k,d,alpha.d)
            exp( - (lam1-lam2)*alpha.d * x^d + k * (log(lam1) - log(lam2)))
        }
        NNdist <- function(X, k)
        {
            ##--------------------------------------------------------------
            ## nearest-neighbor via S's  dist()
            ##
            n <- nrow(X)
            distances <- dist(X)

            ##  This next part sorts through the R distance object and forms kNNd,
            ##  kth nearest neighbour distance, for each point.
            ##
            kNNd <- numeric(n)
            N <- (n - 1):0
            I <- c(0, cumsum(N[-1]))
            J <- c(0, I + n - 1)
            a <- z <- NULL
            for(j in 1:n) {
                if(j > 1)
                    a <- i + I[1:i]
                if(j < n)
                    z <- J[j] + 1:N[j]
                kNNd[j] <- sort(distances[c(a, z)], partial = k)[k]
                i <- j
            }
            kNNd
        }
        ## end of NNdist

        ##--  begin{nclean.sub} -------------------------------------
        stopifnot(length(d <- dim(X)) == 2, d >= 1)
        n <- d[1]
        d <- d[2]
        ##
        alpha.d <- (2 * pi^(d/2))/(d * gamma(d/2))
        ##
        ## Now use k.dists in E-M algorithm, first get starting guesses.
        ##
        delta <- integer(n)
        delta[k.dists > (min(k.dists) + diff(range(k.dists))/3)] <- 1L
        p <- 0.5
        lambda1 <- k/(alpha.d * mean(k.dists[delta == 0L]^d))
        lambda2 <- k/(alpha.d * mean(k.dists[delta == 1L]^d))
        loglik.old <- 0
        loglik.new <- 1
        ##
        ## Iterator starts here,
        ##
        it <- 0
        while(abs(loglik.new - loglik.old)/(1+abs(loglik.new)) > convergence)
        {
            it <- it+1
            ## E - step
            if(FALSE) { ## old code (suffering from "0/0")
                dDk.l1 <- dDk(k.dists, lambda1, k = k, d = d, alpha.d = alpha.d)
                dDk.l2 <- dDk(k.dists, lambda2, k = k, d = d, alpha.d = alpha.d)
                delta <- (p * dDk.l1) / (p* dDk.l1 + (1 - p) * dDk.l2)
                delta0 <-(1) / (1 + (1/p - 1) * dDk.l2/dDk.l1)
            }
            delta <- (1) / (1 + (1/p - 1) * q.dDk(k.dists, lambda2, lambda1,
                                                  k = k, d = d, alpha.d = alpha.d))
            ## M - step
            p <- sum(delta)/n
            akNN.d <- alpha.d * k.dists^d
            lambda1 <- k * sum(delta)      / sum(a..d <- akNN.d * delta)
            lambda2 <- k * sum((1 - delta))/ sum(a.1d <- akNN.d * (1 - delta))
            loglik.old <- loglik.new
            loglik.new <- sum(   - p    * lambda1 * a..d
                              - (1 - p) * lambda2 * a.1d +
                                 delta    * k * log(lambda1 * alpha.d) +
                              (1 - delta) * k * log(lambda2 * alpha.d))
        }
        ##
        ## z will be the classifications. 1= in cluster. 0= in noise.
        ##

        probs <- 1 / (1 + q.dDk(k.dists, lambda2, lambda1,
                                k = k, d = d, alpha.d = alpha.d))
        mprob <- 1. - probs
        mu1 <- colSums(probs * X)/sum(probs)
        mu2 <- colSums(mprob * X)/sum(mprob)
        tpsig1 <- t(X) - mu1
        tpsig2 <- t(X) - mu2
        Sig1 <- tpsig1 %*% (probs * t(tpsig1))/sum(probs)
        Sig2 <- tpsig2 %*% (mprob * t(tpsig2))/sum(mprob)
        sd1 <- sqrt(diag(Sig1))
        sd2 <- sqrt(diag(Sig2))
        zz <- list(z = round(probs), kthNND = k.dists, probs = probs,
                   p = p, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2,
                   lambda1 = lambda1, lambda2 = lambda2, Sig1 = Sig1, Sig2 = Sig2,
                   iter = it)
        return(zz)
    }
    ##----- end of nclean.sub-----------------------------------------------

    ##---begin{cov.nnve} ---{main}-----------------------------------------
    ##
    X <- as.matrix(X)
    n <- nrow(X)
    d <- ncol(X)
    stopifnot(n >= 3, d >= 1)

    S.mean <- apply(X, 2, median)
    S.sd <- apply(X, 2, mad)
    my.scale <- function(x) {
        sweep(sweep(x, 2, S.mean, check.margin = FALSE),
              2, S.sd, "/", check.margin = FALSE)
    }
    X. <- my.scale(X)
    ##
    ##  NNC based on original data
    ##
    orgNNC <- nclean.sub(X, k, convergence = emconv, s.X = X.)
    nnoise <- min(c(sum(1 - orgNNC$z), round(pnoise * n)))
    knnd <- orgNNC$kthNND
    ord <- (n + 1) - rank(knnd)
    muT <- orgNNC$mu1
    SigT <- orgNNC$Sig1
    SigT <- (SigT + t(SigT))/2.
    SigTN <- diag(orgNNC$sd1^2, d)
    if(nnoise > 6) {
        ncho <- nnoise
        ncho1 <- floor(ncho/2)
        ncho2 <- ncho - ncho1
        cho <- (1:n)[ord <= ncho1]
        xcho <- X[cho, , drop=FALSE]
        ev <- eigen(SigT)
        evv <- ev$values
        minv <- max((1:d)[evv > 1e-12])
        if(minv > 2) {
            vv1 <- ev$vectors[, (minv - 1)]
            vv2 <- ev$vectors[, minv]
        }
        else {
            vv1 <- ev$vectors[, 1]
            vv2 <- ev$vectors[, 2]
        }
        ot <- acos(sum(vv1 * vv2)/(sum(vv1^2) * sum(vv2^2))^0.5)
        for(kk1 in 1:(ncho2)) {
            pseg <- kk1/(ncho2 + 1) * ot
            xcho <- rbind(xcho, (sin(pseg) * vv1 +
                                 cos(pseg) * vv2 + muT))
        }
    }
    else {
        nnoise <- 3
        cho <- (1:n)[ord <= nnoise]
        xcho <- X[cho, , drop=FALSE]
    }
    n2 <- nrow(xcho)
    schox <- mahalanobis(xcho, muT, SigTN)
    Nc <- matrix(rep(muT, each=n2), nrow = n2)
    Ndir <- (xcho - Nc)/(schox^0.5)
    ##
    ## initial set up
    ##
    ch1 <- qchisq(c(.01, .0001), d, lower.tail=FALSE)
    Xa <- seq(ch1[1], ch1[2], length = 6)
    gap <- Xa[2] - Xa[1]
    initv <- diag(orgNNC$Sig1)
    OldP <- orgNNC$probs
    SaveL <- list(orgNNC[c("mu1","Sig1")])
    SaveP <- list(OldP)
    ## __FIXME__ From here on + while() { .. } below, have repeat {...}
    xa <- Xa[1]
    Np <- Nc - Ndir * (xa^0.5)
    X.N <- rbind(X, Np)
    sXN <- rbind(X., my.scale(Np))
    updNNC <- nclean.sub(X.N, k, convergence = emconv, s.X = sXN)
    SaveL <- c(SaveL, list(updNNC[c("mu1","Sig1")]))
    SaveP <- c(SaveP, list(updNNC$probs[1:n])) # have length (n + n2)
    ##sda<-Mtovec(orgNNC$Sig1)
    ##sda save the results corresponding to xa = qchisq(.99,d)
    stopv <- diag(updNNC$Sig1)
    time1 <- 2
    ##
    ##
    while((time1 <= 6) && all(stopv < (1 + bound) * initv)) {
        xa <- Xa[time1]
        Np <- Nc - Ndir * (xa^0.5)
        X.N[n+ 1:n2,] <- Np
        sXN[n+ 1:n2,] <- my.scale(Np)
        updNNC <- nclean.sub(X.N, k, convergence = emconv, s.X = sXN)
        SaveL <- c(SaveL, list(updNNC[c("mu1","Sig1")]))
        SaveP <- c(SaveP[2], list(updNNC$probs[1:n])) # always keep the last two
        time1 <- time1 + 1
        stopv <- diag(updNNC$Sig1)
        NULL
    }
    ##
    ##   Procedure stop if the added noise cause a "surge" within the range
    ##
    nS <- length(SaveL)
    if(all(stopv < (1 + bound) * initv)) {
        ans <- SaveL[[nS]]
        NewP <- SaveP[[2]]
        ##
        ##  adding extension
        ##
        if(extension) {
            Fstop <- FALSE
            tpv <- stopv
            while(!Fstop && all(stopv < (1 + bound) * initv)) {
                xa <- xa + gap
                startv <- stopv
                Np <- Nc - Ndir * (xa^0.5)
                X.N[n+ 1:n2,] <- Np
                sXN[n+ 1:n2,] <- my.scale(Np)
                updNNC <- nclean.sub(X.N, k, convergence = emconv, s.X = sXN)
                SaveL <- c(SaveL, list(updNNC[c("mu1","Sig1")]))
                SaveP <- c(SaveP[2], list(updNNC$probs[1:n])) # always keep the last two
                ntpv <- diag(updNNC$Sig1)
                stopv <- colMeans(rbind(startv * 2 - tpv, ntpv))
                tpv <- ntpv
                Fstop <- all(abs(stopv - startv) <= (1+abs(startv)) * devsm)
            } ## end-while

            nS <- length(SaveL)
            ## checking the stop criterion at the end of extension
            if(all(stopv < (1 + bound) * initv)) {
                ans <- SaveL[[nS]]
                NewP <- SaveP[[2]]
            }
            else {
                ans <- SaveL[[nS - 1]]
                NewP <- SaveP[[1]]
            }
        }
    }
    else {
        ans <- SaveL[[nS - 1]]
        NewP <- SaveP[[1]]
    }
    ##
    list(cov = ans[["Sig1"]], mu = ans[["mu1"]],
         postprob = NewP, classification = round(NewP), iter = updNNC$iter,
         innc =
         list(cov = orgNNC$Sig1, mu = orgNNC$mu1, iter = orgNNC$iter,
              postprob = OldP, classification = round(OldP)))
}

cov.nnve <- function(datamat, k = 12, ...)
{
    .Deprecated("covNNC")
    covNNC(datamat, k=k, ...)
}
