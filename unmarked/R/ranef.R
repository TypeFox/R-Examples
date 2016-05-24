



# ----------------- Empirical Bayes Methods ------------------------------



setGeneric("ranef",
    function(object, ...) standardGeneric("ranef"))


setClass("unmarkedRanef",
    representation(post = "array"))


setClassUnion("unmarkedFitGMMorGDS",
              c("unmarkedFitGMM", "unmarkedFitGDS"))


setMethod("ranef", "unmarkedFitPCount",
    function(object, ...)
{
    lam <- predict(object, type="state")[,1] # Too slow
    R <- length(lam)
    p <- getP(object)
    K <- object@K
    N <- 0:K
    y <- getY(getData(object))
    srm <- object@sitesRemoved
    if(length(srm) > 0)
        y <- y[-object@sitesRemoved,]
    post <- array(NA_real_, c(R, length(N), 1))
    colnames(post) <- N
    mix <- object@mixture
    for(i in 1:R) {
        switch(mix,
               P  = f <- dpois(N, lam[i], log=TRUE),
               NB = {
                   alpha <- exp(coef(object, type="alpha"))
                   f <- dnbinom(N, mu=lam[i], size=alpha, log=TRUE)
               },
               ZIP = {
                   psi <- plogis(coef(object, type="psi"))
                   f <- (1-psi)*dpois(N, lam[i])
                   f[1] <- psi + (1-psi)*exp(-lam[i])
                   f <- log(f)
               })
        g <- rep(0, K+1)
        for(j in 1:ncol(y)) {
            if(is.na(y[i,j]) | is.na(p[i,j]))
                next
            g <- g + dbinom(y[i,j], N, p[i,j], log=TRUE)
        }
        fudge <- exp(f+g)
        post[i,,1] <- fudge / sum(fudge)
    }
    new("unmarkedRanef", post=post)
})






setMethod("ranef", "unmarkedFitOccu",
    function(object, ...)
{
    psi <- predict(object, type="state")[,1]
    R <- length(psi)
    p <- getP(object)
    z <- 0:1
    y <- getY(getData(object))
    srm <- object@sitesRemoved
    if(length(srm) > 0)
        y <- y[-object@sitesRemoved,]
    y[y>1] <- 1
    post <- array(0, c(R,2,1))
    colnames(post) <- z
    for(i in 1:R) {
        f <- dbinom(z, 1, psi[i])
        g <- rep(1, 2)
        for(j in 1:ncol(y)) {
            if(is.na(y[i,j]) | is.na(p[i,j]))
                next
            g <- g * dbinom(y[i,j], 1, z*p[i,j])
        }
        fudge <- f*g
        post[i,,1] <- fudge / sum(fudge)
    }
    new("unmarkedRanef", post=post)
})


setMethod("ranef", "unmarkedFitOccuFP", function(object, na.rm = FALSE)
{
  cat("ranef is not implemented for occuFP at this time")
})






setMethod("ranef", "unmarkedFitOccuRN",
    function(object, K, ...)
{
    if(missing(K)) {
        warning("You did not specify K, the maximum value of N, so it was set to 50")
        K <- 50
    }
    lam <- predict(object, type="state")[,1] # Too slow
    R <- length(lam)
    r <- getP(object)
    N <- 0:K
    y <- getY(getData(object))
    srm <- object@sitesRemoved
    if(length(srm) > 0)
        y <- y[-object@sitesRemoved,]
    y[y>1] <- 1
    post <- array(NA_real_, c(R, length(N), 1))
    colnames(post) <- N
    for(i in 1:R) {
        f <- dpois(N, lam[i])
        g <- rep(1, K+1)
        for(j in 1:ncol(y)) {
            if(is.na(y[i,j]) | is.na(r[i,j]))
                next
            p.ijn <- 1 - (1-r[i,j])^N
            g <- g * dbinom(y[i,j], 1, p.ijn)
        }
        fudge <- f*g
        post[i,,1] <- fudge / sum(fudge)
    }
    new("unmarkedRanef", post=post)
})





setMethod("ranef", "unmarkedFitMPois",
    function(object, K, ...)
{
    y <- getY(getData(object))
    srm <- object@sitesRemoved
    if(length(srm) > 0)
        y <- y[-object@sitesRemoved,]
    if(missing(K)) {
        warning("You did not specify K, the maximum value of N, so it was set to max(y)+50")
        K <- max(y, na.rm=TRUE)+50
    }

    lam <- predict(object, type="state")[,1]
    R <- length(lam)
    cp <- getP(object)
    cp <- cbind(cp, 1-rowSums(cp))
    N <- 0:K
    post <- array(0, c(R, K+1, 1))
    colnames(post) <- N
    for(i in 1:R) {
        f <- dpois(N, lam[i])
        g <- rep(1, K+1)
        if(any(is.na(y[i,])) | any(is.na(cp[i,])))
            next
        for(k in 1:(K+1)) {
            yi <- y[i,]
            ydot <- N[k] - sum(yi)
            if(ydot<0) {
                g[k] <- 0
                next
            }
            yi <- c(yi, ydot)
            g[k] <- g[k] * dmultinom(yi, size=N[k], prob=cp[i,])
        }
        fudge <- f*g
        post[i,,1] <- fudge / sum(fudge)
    }
    new("unmarkedRanef", post=post)
})







setMethod("ranef", "unmarkedFitDS",
    function(object, K, ...)
{
    y <- getY(getData(object))
    srm <- object@sitesRemoved
    if(length(srm) > 0)
        y <- y[-object@sitesRemoved,]
    if(missing(K)) {
        warning("You did not specify K, the maximum value of N, so it was set to max(y)+50")
        K <- max(y, na.rm=TRUE)+50
    }
    lam <- predict(object, type="state")[,1]
    R <- length(lam)
    J <- ncol(y)
    survey <- object@data@survey
    tlength <- object@data@tlength
    db <- object@data@dist.breaks
    w <- diff(db)
    unitsIn <- object@data@unitsIn
    unitsOut <- object@unitsOut
    if(identical(object@output, "density")) {
        a <- matrix(NA, R, J)
        switch(survey, line = {
            for (i in 1:R) {
                a[i, ] <- tlength[i] * w
            }
        }, point = {
            for (i in 1:R) {
                a[i, 1] <- pi * db[2]^2
                for (j in 2:J)
                    a[i, j] <- pi * db[j + 1]^2 - sum(a[i, 1:(j - 1)])
            }
        })
        switch(survey, line = A <- rowSums(a) * 2, point = A <- rowSums(a))
        switch(unitsIn, m = A <- A/1e+06, km = A <- A)
        switch(unitsOut, ha = A <- A * 100, kmsq = A <- A)
        lam <- lam*A
    }
    cp <- getP(object)
    cp <- cbind(cp, 1-rowSums(cp))
    N <- 0:K
    post <- array(0, c(R, K+1, 1))
    colnames(post) <- N
    for(i in 1:R) {
        f <- dpois(N, lam[i])
        g <- rep(1, K+1)
        if(any(is.na(y[i,])) | any(is.na(cp[i,])))
            next
        for(k in 1:(K+1)) {
            yi <- y[i,]
            ydot <- N[k] - sum(yi)
            if(ydot<0) {
                g[k] <- 0
                next
            }
            yi <- c(yi, ydot)
            g[k] <- g[k] * dmultinom(yi, size=N[k], prob=cp[i,])
        }
        fudge <- f*g
        post[i,,1] <- fudge / sum(fudge)
    }
    new("unmarkedRanef", post=post)
})







setMethod("ranef", "unmarkedFitGMMorGDS",
    function(object, ...)
{

    data <- object@data
    D <- getDesign(data, object@formula)

    Xlam <- D$Xlam
    Xphi <- D$Xphi
    Xdet <- D$Xdet
    y <- D$y  # MxJT
    Xlam.offset <- D$Xlam.offset
    Xphi.offset <- D$Xphi.offset
    Xdet.offset <- D$Xdet.offset
    if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
    if(is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))
    if(is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))
    beta.lam <- coef(object, type="lambda")
    beta.phi <- coef(object, type="phi")
    beta.det <- coef(object, type="det")

    lambda <- exp(Xlam %*% beta.lam + Xlam.offset)
    if(is.null(beta.phi))
        phi <- rep(1, nrow(Xphi))
    else
        phi <- plogis(Xphi %*% beta.phi + Xphi.offset)

    cp <- getP(object)
    cp[is.na(y)] <- NA

    K <- object@K
    M <- 0:K

    nSites <- nrow(y)
    T <- data@numPrimary
    R <- numY(data) / T
    J <- obsNum(data) / T

    phi <- matrix(phi, nSites, byrow=TRUE)

    if(identical(class(object)[1], "unmarkedFitGDS")) {
        if(identical(object@output, "density")) {
            survey <- object@data@survey
            tlength <- object@data@tlength
            unitsIn <- object@data@unitsIn
            unitsOut <- object@unitsOut
            db <- object@data@dist.breaks
            w <- diff(db)
            u <- a <- matrix(NA, nSites, R)
            switch(survey,
                   line = {
                       for(i in 1:nSites) {
                           a[i,] <- tlength[i] * w
                           u[i,] <- a[i,] / sum(a[i,])
                       }
                   },
                   point = {
                       for(i in 1:nSites) {
                           a[i, 1] <- pi*db[2]^2
                           for(j in 2:R)
                               a[i, j] <- pi*db[j+1]^2 - sum(a[i, 1:(j-1)])
                           u[i,] <- a[i,] / sum(a[i,])
                       }
                   })
            switch(survey,
                   line = A <- rowSums(a) * 2,
                   point = A <- rowSums(a))
            switch(unitsIn,
                   m = A <- A / 1e6,
                   km = A <- A)
            switch(unitsOut,
                   ha = A <- A * 100,
                   kmsq = A <- A)
            lambda <- lambda*A # Density to abundance
        }
    }

    cpa <- array(cp, c(nSites,R,T))
    ya <- array(y, c(nSites,R,T))
#    ym <- apply(ya, c(1,3), sum, na.rm=TRUE)

    post <- array(0, c(nSites, K+1, 1))
    colnames(post) <- M
    mix <- object@mixture
    if(identical(mix, "NB"))
        alpha <- exp(coef(object, type="alpha"))
    for(i in 1:nSites) {
        switch(mix,
               P  = f <- dpois(M, lambda[i]),
               # FIXME: Add ZIP
               NB = f <- dnbinom(M, mu=lambda[i], size=alpha))
        g <- rep(1, K+1) # outside t loop
        for(t in 1:T) {
            if(all(is.na(ya[i,,t])) | is.na(phi[i,t]))
                next
            for(k in 1:(K+1)) {
                y.it <- ya[i,,t]
                ydot <- M[k]-sum(y.it, na.rm=TRUE)
                y.it <- c(y.it, ydot)
                if(ydot < 0) {
                    g[k] <- 0
                    next
                }
                cp.it <- cpa[i,,t]*phi[i,t]
                cp.it <- c(cp.it, 1-sum(cp.it, na.rm=TRUE))
                na.it <- is.na(cp.it)
                y.it[na.it] <- NA
                g[k] <- g[k]*dmultinom(y.it[!na.it], M[k], cp.it[!na.it])
            }
        }
        fudge <- f*g
        post[i,,1] <- fudge/sum(fudge)
    }
    new("unmarkedRanef", post=post)
})






setMethod("ranef", "unmarkedFitGPC",
    function(object, ...)
{
    data <- object@data
    D <- getDesign(data, object@formula)

    Xlam <- D$Xlam
    Xphi <- D$Xphi
    Xdet <- D$Xdet
    y <- D$y  # MxJT
    Xlam.offset <- D$Xlam.offset
    Xphi.offset <- D$Xphi.offset
    Xdet.offset <- D$Xdet.offset
    if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
    if(is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))
    if(is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))

    beta.lam <- coef(object, type="lambda")
    beta.phi <- coef(object, type="phi")
    beta.det <- coef(object, type="det")

    R <- nrow(y)
    T <- data@numPrimary
    J <- ncol(y) / T

    lambda <- exp(Xlam %*% beta.lam + Xlam.offset)
    if(is.null(beta.phi))
        phi <- rep(1, nrow(Xphi))
    else
        phi <- plogis(Xphi %*% beta.phi + Xphi.offset)
    phi <- matrix(phi, R, byrow=TRUE)

    p <- getP(object)
    p[is.na(y)] <- NA
    pa <- array(p, c(R,J,T))
    ya <- array(y, c(R,J,T))

    K <- object@K
    M <- N <- 0:K
    lM <- K+1

    post <- array(0, c(R, K+1, 1))
    colnames(post) <- M
    mix <- object@mixture
    if(identical(mix, "NB"))
        alpha <- exp(coef(object, type="alpha"))
    for(i in 1:R) {
        switch(mix,
               P  = f <- dpois(M, lambda[i]),
               # FIXME: Add ZIP
               NB = f <- dnbinom(M, mu=lambda[i], size=alpha))
        ghi <- rep(0, lM)
        for(t in 1:T) {
            gh <- matrix(-Inf, lM, lM)
            for(m in M) {
                if(m < max(ya[i,,], na.rm=TRUE)) {
                    gh[,m+1] <- -Inf
                    next
                }
                if(is.na(phi[i,t])) {
                    g <- rep(0, lM)
                    g[N>m] <- -Inf
                }
                else
                    g <- dbinom(N, m, phi[i,t], log=TRUE)
                h <- rep(0, lM)
                for(j in 1:J) {
                    if(is.na(ya[i,j,t]) | is.na(pa[i,j,t]))
                        next
                    h <- h + dbinom(ya[i,j,t], N, pa[i,j,t], log=TRUE)
                }
                gh[,m+1] <- g + h
            }
            ghi <- ghi + log(colSums(exp(gh)))
        }
        fgh <- exp(f + ghi)
        prM <- fgh/sum(fgh)
        post[i,,1] <- prM
    }
    new("unmarkedRanef", post=post)
})






setMethod("ranef", "unmarkedFitColExt",
    function(object, ...)
{

#    stop("method not written")

    data <- object@data
    M <- numSites(data)
    nY <- data@numPrimary
    J <- obsNum(data)/nY
    psiParms <- coef(object, 'psi')
    detParms <- coef(object, 'det')
    colParms <- coef(object, 'col')
    extParms <- coef(object, 'ext')
    formulaList <- list(psiformula=object@psiformula,
        gammaformula=object@gamformula,
        epsilonformula=object@epsformula,
        pformula=object@detformula)
    designMats <- getDesign(object@data, object@formula)
    V.itj <- designMats$V
    X.it.gam <- designMats$X.gam
    X.it.eps <- designMats$X.eps
    W.i <- designMats$W

#    yumf <- getY(object@data)
#    yumfa <- array(yumf, c(M, J, nY))

    y <- designMats$y
    y[y>1] <- 1
    ya <- array(y, c(M, J, nY))

    psiP <- plogis(W.i %*% psiParms)
    detP <- plogis(V.itj %*% detParms)
    colP <- plogis(X.it.gam  %*% colParms)
    extP <- plogis(X.it.eps %*% extParms)

    detP <- array(detP, c(J, nY, M))
    colP <- matrix(colP, M, nY, byrow = TRUE)
    extP <- matrix(extP, M, nY, byrow = TRUE)

    ## create transition matrices (phi^T)
    phis <- array(NA,c(2,2,nY-1,M)) #array of phis for each
    for(i in 1:M) {
        for(t in 1:(nY-1)) {
            phis[,,t,i] <- matrix(c(1-colP[i,t], colP[i,t], extP[i,t],
                1-extP[i,t]))
            }
        }

    ## first compute latent probs
    x <- array(NA, c(2, nY, M))
    x[1,1,] <- 1-psiP
    x[2,1,] <- psiP
    for(i in 1:M) {
        for(t in 2:nY) {
            x[,t,i] <- (phis[,,t-1,i] %*% x[,t-1,i])
            }
        }

    z <- 0:1
    post <- array(NA_real_, c(M, 2, nY))
    colnames(post) <- z

    for(i in 1:M) {
        for(t in 1:nY) {
            g <- rep(1, 2)
            for(j in 1:J) {
                if(is.na(ya[i,j,t]) | is.na(detP[j,t,i]))
                    next
                g <- g * dbinom(ya[i,j,t], 1, z*detP[j,t,i])
            }
            tmp <- x[,t,i] * g
            post[i,,t] <- tmp/sum(tmp)
        }
    }

    new("unmarkedRanef", post=post)
})







setMethod("ranef", "unmarkedFitPCO",
    function(object, ...)
{
#    browser()
    dyn <- object@dynamics
    formlist <- object@formlist
    formula <- as.formula(paste(unlist(formlist), collapse=" "))
    D <- getDesign(object@data, formula)
    delta <- D$delta
    deltamax <- max(delta, na.rm=TRUE)
    if(!.hasSlot(object, "immigration")) #For backwards compatibility
      imm <- FALSE
    else
      imm <- object@immigration

    lam <- predict(object, type="lambda")[,1] # Slow, use D$Xlam instead
    R <- length(lam)
    T <- object@data@numPrimary
    p <- getP(object)
    K <- object@K
    N <- 0:K
    y <- getY(getData(object))
    J <- ncol(y)/T
    if(dyn != "notrend") {
        gam <- predict(object, type="gamma")[,1]
        gam <- matrix(gam, R, T-1, byrow=TRUE)
    }
    if(!identical(dyn, "trend")) {
        om <- predict(object, type="omega")[,1]
        om <- matrix(om, R, T-1, byrow=TRUE)
    }
    else
        om <- matrix(0, R, T-1)
    if(imm) {
        iota <- predict(object, type="iota")[,1]
        iota <- matrix(iota, R, T-1, byrow=TRUE)
    }
    else
      iota <- matrix(0, R, T-1)
    srm <- object@sitesRemoved
    if(length(srm) > 0)
        y <- y[-object@sitesRemoved,]
    ya <- array(y, c(R, J, T))
    pa <- array(p, c(R, J, T))
    post <- array(NA_real_, c(R, length(N), T))
    colnames(post) <- N
    mix <- object@mixture
    if(dyn=="notrend")
        gam <- lam*(1-om)

    if(dyn %in% c("constant", "notrend")) {
        tp <- function(N0, N1, gam, om, iota) {
            c <- 0:min(N0, N1)
            sum(dbinom(c, N0, om) * dpois(N1-c, gam))
        }
    } else if(dyn=="autoreg") {
        tp <- function(N0, N1, gam, om, iota) {
            c <- 0:min(N0, N1)
            sum(dbinom(c, N0, om) * dpois(N1-c, gam*N0 + iota))
        }
    } else if(dyn=="trend") {
        tp <- function(N0, N1, gam, om, iota) {
            dpois(N1, gam*N0 + iota)
        }
    } else if(dyn=="ricker") {
        tp <- function(N0, N1, gam, om, iota) {
            dpois(N1, N0*exp(gam*(1-N0/om)) + iota)
        }
    } else if(dyn=="gompertz") {
        tp <- function(N0, N1, gam, om, iota) {
            dpois(N1, N0*exp(gam*(1-log(N0 + 1)/log(om + 1))) + iota)
        }
    }
    for(i in 1:R) {
        P <- matrix(1, K+1, K+1)
        switch(mix,
               P  = g2 <- dpois(N, lam[i]),
               NB = {
                   alpha <- exp(coef(object, type="alpha"))
                   g2 <- dnbinom(N, mu=lam[i], size=alpha)
               },
               ZIP = {
                   psi <- plogis(coef(object, type="psi"))
                   g2 <- (1-psi)*dpois(N, lam[i])
                   g2[1] <- psi + (1-psi)*exp(-lam[i])
               })
        g1 <- rep(1, K+1)
        for(j in 1:J) {
            if(is.na(ya[i,j,1]) | is.na(pa[i,j,1]))
                next
            g1 <- g1 * dbinom(ya[i,j,1], N, pa[i,j,1])
        }
        g1g2 <- g1*g2
        post[i,,1] <- g1g2 / sum(g1g2)
        for(t in 2:T) {
            if(!is.na(gam[i,t-1]) & !is.na(om[i,t-1])) {
                for(n0 in N) {
                    for(n1 in N) {
                        P[n0+1, n1+1] <- tp(n0, n1, gam[i,t-1], om[i,t-1], iota[i,t-1])
                    }
                }
            }
            delta.it <- delta[i,t-1]
            if(delta.it > 1) {
                P1 <- P
                for(d in 2:delta.it) {
                    P <- P %*% P1
                }
            }
            g1 <- rep(1, K+1)
            for(j in 1:J) {
                if(is.na(ya[i,j,t]) | is.na(pa[i,j,t]))
                    next
                g1 <- g1 * dbinom(ya[i,j,t], N, pa[i,j,t])
            }
            g <- colSums(P * post[i,,t-1]) * g1
            post[i,,t] <- g / sum(g)
        }
    }
    new("unmarkedRanef", post=post)
})






setGeneric("bup", function(object, stat=c("mean", "mode"), ...)
    standardGeneric("bup"))
setMethod("bup", "unmarkedRanef",
          function(object, stat=c("mean", "mode"), ...) {
    stat <- match.arg(stat)
    post <- object@post
    re <- as.integer(colnames(post))
    if(identical(stat, "mean"))
        out <- apply(post, c(1,3), function(x) sum(re*x))
    else if(identical(stat, "mode"))
        out <- apply(post, c(1,3), function(x) re[which.max(x)])
    out <- drop(out)
    return(out)
})








setMethod("confint", "unmarkedRanef", function(object, parm, level=0.95)
{
    if(!missing(parm))
        warning("parm argument is ignored. Did you mean to specify level?")
    post <- object@post
    N <- as.integer(colnames(post))
    R <- nrow(post)
    T <- dim(post)[3]
    CI <- array(NA_real_, c(R,2,T))
    alpha <- 1-level
    c1 <- alpha/2
    c2 <- 1-c1
    colnames(CI) <- paste(c(c1,c2)*100, "%", sep="")
    for(i in 1:R) {
        for(t in 1:T) {
            pr <- post[i,,t]
            ed <- cumsum(pr)
            lower <- N[which(ed >= c1)][1]
            upper <- N[which(ed >= c2)][1]
            CI[i,,t] <- c(lower, upper)
        }
    }
    CI <- drop(CI) # Convert to matrix if T==1
    return(CI)
})






setMethod("show", "unmarkedRanef", function(object)
{
    post <- object@post
    dims <- dim(post)
    T <- dims[3]
    if(T==1)
        print(cbind(Mean=bup(object, stat="mean"),
                    Mode=bup(object, stat="mode"), confint(object)))
    else if(T>1) {
        means <- bup(object, stat="mean")
        modes <- bup(object, stat="mode")
        CI <- confint(object)
        out <- array(NA_real_, c(dims[1], 4, T))
        dimnames(out) <- list(NULL,
                              c("Mean",
                                "Mode", "2.5%", "97.5%"),
                              paste("Year", 1:T, sep=""))
        for(t in 1:T) {
            out[,,t] <- cbind(means[,t],
                              modes[,t], CI[,,t])
        }
        print(out)
    }
})



setAs("unmarkedRanef", "array", function(from) {
    post <- from@post
    dims <- dim(post)
    R <- dims[1]
    T <- dims[3]
    dimnames(post) <- list(1:R, colnames(post), 1:T)
    post <- drop(post)
    return(post)
})


setAs("unmarkedRanef", "data.frame", function(from) {
    post <- from@post
    dims <- dim(post)
    R <- dims[1]
    lN <- dims[2]
    T <- dims[3]
    N <- as.integer(colnames(post))
    N.ikt <- rep(rep(N, each=R), times=T)
    site <- rep(1:R, times=lN*T)
    year <- rep(1:T, each=R*lN)
    dat <- data.frame(site=site, year=year, N=N.ikt, p=as.vector(post))
    dat <- dat[order(dat$site),]
    if(T==1)
        dat$year <- NULL
    return(dat)
})



setMethod("plot", c("unmarkedRanef", "missing"), function(x, y, ...)
{
    post <- x@post
    T <- dim(post)[3]
    N <- as.integer(colnames(post))
    xlb <- ifelse(length(N)>2, "Abundance", "Occurrence")
    ylb <- "Probability"
    dat <- as(x, "data.frame")
    site.c <- as.character(dat$site)
    nc <- nchar(site.c)
    mc <- max(nc)
    dat$site.c <- paste("site", sapply(site.c, function(x)
         paste(paste(rep("0", mc-nchar(x)), collapse=""), x, sep="")),
         sep="")
    if(T==1)
        xyplot(p ~ N | site.c, dat, type="h", xlab=xlb, ylab=ylb, ...)
    else if(T>1) {
        year.c <- as.character(dat$year)
        nc <- nchar(year.c)
        mc <- max(nc)
        dat$year.c <- paste("year", sapply(year.c, function(x)
            paste(paste(rep("0", mc-nchar(x)), collapse=""), x, sep="")),
            sep="")
        xyplot(p ~ N | site.c+year.c, dat, type="h", xlab=xlb, ylab=ylb, ...)
    }
})












