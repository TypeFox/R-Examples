## TO DO
##
##   - 'best', 'exact' and 'deterministic' options for nsamp, as in CovMcd
##
CovSest <- function(x,
                    bdp=0.5,
                    arp=0.1,
                    eps=1e-5,
                    maxiter=120,
                    nsamp=500,
                    seed=NULL,
                    trace=FALSE,
                    tolSolve=1e-14,
                    scalefn,
                    maxisteps=200,
                    initHsets = NULL, save.hsets = FALSE,
                    method=c("sfast", "surreal", "bisquare", "rocke", "suser", "sdet"),
                    control,
                    t0,
                    S0,
                    initcontrol
                )
{

    ## Analize and validate the input parameters ...
    ## if a control object was supplied, take the option parameters from it,
    ## but if single parameters were passed (not defaults) they will override the
    ## control object.

    method <- match.arg(method)
    if(!missing(control)){
        defcontrol <- CovControlSest()     # default control
        if(bdp == defcontrol@bdp)           bdp <- control@bdp          # for s-fast and surreal
        if(arp == defcontrol@arp)           arp <- control@arp          # for rocke type
        if(eps == defcontrol@eps)           eps <- control@eps          # for bisquare and rocke
        if(maxiter == defcontrol@maxiter)   maxiter <- control@maxiter  # for bisquare and rocke

        if(nsamp == defcontrol@nsamp)       nsamp <- control@nsamp
        if(is.null(seed) || seed == defcontrol@seed)         seed <- control@seed
        if(trace == defcontrol@trace)       trace <- control@trace
        if(tolSolve == defcontrol@tolSolve) tolSolve <- control@tolSolve
        if(method == defcontrol@method)     method <- control@method
    }

    if(length(seed) > 0) {
        if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))  {
            seed.keep <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir=.GlobalEnv))
        }
        assign(".Random.seed", seed, envir=.GlobalEnv)
    }

    if(is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
                    dimnames = list(names(x), deparse(substitute(x))))

    xcall <- match.call()

    ## drop all rows with missing values (!!) :
    na.x <- !is.finite(x %*% rep(1, ncol(x)))
    ok <- !na.x
    x <- x[ok, , drop = FALSE]
    dx <- dim(x)
    if(!length(dx))
        stop("All observations have missing values!")
    dimn <- dimnames(x)
    n <- dx[1]
    p <- dx[2]

    if(n <= p + 1)
        stop(if (n <= p) "n <= p -- you can't be serious!" else "n == p+1  is too small sample size")
    if(n < 2 * p)
    { ## p+1 < n < 2p
        warning("n < 2 * p, i.e., possibly too small sample size")
    }

    if(method == "surreal" && nsamp == 500)     # default
        nsamp = 600*p

    ## compute the constants c1 and kp
    ## c1 = Tbsc(bdp, p)
    ## kp = (c1/6) * Tbsb(c1, p)
    cobj <- .csolve.bw.S(bdp, p)
    cc <- cobj$cc
    kp <- cobj$kp

    if(trace)
        cat("\nFAST-S...: bdp, p, cc, kp=", bdp, p, cc, kp, "\n")

    if(missing(scalefn) || is.null(scalefn))
    {
        scalefn <- if(n <= 1000) Qn else scaleTau2
    }

    mm <- if(method == "sfast")         ..CSloc(x, nsamp=nsamp, kp=kp, cc=cc, trace=trace)
          else if(method == "suser")    ..fastSloc(x, nsamp=nsamp, kp=kp, cc=cc, trace=trace)
          else if(method == "surreal")  ..covSURREAL(x, nsamp=nsamp, kp=kp, c1=cc, trace=trace, tol.inv=tolSolve)
          else if(method == "bisquare") ..covSBic(x, arp, eps, maxiter, t0, S0, nsamp, seed, initcontrol, trace=trace)
          else if(method == "rocke")    ..covSRocke(x, arp, eps, maxiter, t0, S0, nsamp, seed, initcontrol, trace=trace)
          else                          ..detSloc(x, hsets.init=initHsets, save.hsets=save.hsets, scalefn=scalefn, kp=kp, cc=cc, trace=trace)

    ans <- new("CovSest",
               call = xcall,
               iter=mm$iter,
               crit=mm$crit,
               cov=mm$cov,
               center=mm$center,
               n.obs=n,
               iBest=0,
               cc=mm$cc,
               kp=mm$kp,
               X = as.matrix(x),
               method=mm$method)

    if(method == "sdet")
    {
        ans@iBest <- mm$iBest
        ans@nsteps <- mm$hset.csteps
        if(save.hsets)
            ans@initHsets <- mm$initHsets
    }
    ans
}


##
## Compute the constants kp and c1 for the Tukey Biweight rho-function for S
##
.csolve.bw.S <- function(bdp, p)
{
    ## Compute the constant kp (used in Tbsc)
    Tbsb <- function(c, p)
    {
        ksiint <- function(c, s, p)
        {
            (2^s) * gamma(s + p/2) * pgamma(c^2/2, s + p/2)/gamma(p/2)
        }

        y1 = ksiint(c,1,p)*3/c - ksiint(c,2,p)*3/(c^3) + ksiint(c,3,p)/(c^5)
        y2 = c*(1-pchisq(c^2,p))
        return(y1 + y2)
    }

    ## Compute the tunning constant c1 for the  S-estimator with
    ##  Tukey's biweight function given the breakdown point (bdp)
    ##  and the dimension p
    Tbsc <- function(bdp, p)
    {
        eps = 1e-8
        maxit = 1e3
        cnew <- cold <- sqrt(qchisq(1 - bdp, p))

        ## iterate until the change is less than the tollerance or
        ##  max number of iterations reached
        for(iter in 1:maxit)
        {
            cnew <- Tbsb(cold, p)/bdp
            if(abs(cold - cnew) <= eps)
                break
            cold <- cnew
        }
        return(cnew)
    }

    cc = Tbsc(bdp, p)
    kp = (cc/6) * Tbsb(cc, p)

    return(list(cc=cc, kp=kp))
}

##
##  A fast procedure to compute an S-estimator similar to the one proposed
##  for regression by Salibian-Barrera, M. and Yohai, V.J. (2005),
##  "A fast algorithm for S-regression estimates". This current C implemention
##  is by Valentin Todorov.
##
## Input:
##      x      - a data matrix of size (n,p)
##      nsamp  - number of sub-samples (default=20)
##      k      - number of refining iterations in each subsample (default=2)
##      best.r - number of "best betas" to remember from the subsamples.
##                  These will be later iterated until convergence (default=5)
##      kp, cc  - tunning constants for the  S-estimator with Tukey's biweight
##                function given the breakdown point (bdp) and the dimension p
##
## Output: a list with components
##      center  - robust estimate of location (vector: length)
##      cov     - robust estimate of scatter (matrix: p,p)
##      crit    - value of the objective function (number)
##
..CSloc <- function(x, nsamp=500, kstep=2, best.r=5, kp, cc, trace=FALSE)
{
    dimn <- dimnames(x)
    if(!is.matrix(x))
        x = as.matrix(x)

    n <- nrow(x)
    p <- ncol(x)

    maxit <- 200    # max number of iterations
    rtol <- 1e-7    # convergence tolerance

    b <- .C("sest",
	    as.double(x),
	    as.integer(n),
	    as.integer(p),
	    as.integer(nsamp),
        center = as.double(rep(0,p)),
        cov = as.double(rep(0,p*p)),
        scale = as.double(0),
	    as.double(cc),
	    as.double(kp),
	    as.integer(best.r),                    # number of (refined) to be retained for full iteration
	    as.integer(kstep),                     # number of refining steps for each candidate
	    as.integer(maxit),                     # number of (refined) to be retained for full iteration
	    as.double(rtol),
	    converged = logical(1)
	    )

    scale <- b$scale
    center <- b$center
    cov <- matrix(b$cov, nrow=p)

    if(scale < 0)
	   cat("\nC function sest() exited prematurely!\n")

    ## names(b)[names(b) == "k.max"] <- "k.iter" # maximal #{refinement iter.}
    ## FIXME: get 'res'iduals from C

    if(!is.null(nms <- dimn[[2]])) {
        names(center) <- nms
        dimnames(cov) <- list(nms,nms)
    }

    return(list(
               center=as.vector(center),
               cov=cov,
               crit=scale,
               iter=nsamp,
               kp=kp,
               cc=cc,
               method="S-estimates: S-FAST"))
}

##
##  A fast procedure to compute an S-estimator similar to the one proposed
##  for regression by Salibian-Barrera, M. and Yohai, V.J. (2005),
##  "A fast algorithm for S-regression estimates". This version for
##  multivariate location/scatter is adapted from the one implemented
##  by Kristel Joossens, K.U. Leuven, Belgium
##  and Ella Roelant, Ghent University, Belgium.
##
## Input:
##      x      - a data matrix of size (n,p)
##      nsamp  - number of sub-samples (default=20)
##      k      - number of refining iterations in each subsample (default=2)
##      best.r - number of "best betas" to remember from the subsamples.
##                  These will be later iterated until convergence (default=5)
##      kp, cc  - tunning constants for the  S-estimator with Tukey's biweight
##                function given the breakdown point (bdp) and the dimension p
##
## Output: a list with components
##      center  - robust estimate of location (vector: length)
##      cov     - robust estimate of scatter (matrix: p,p)
##      crit    - value of the objective function (number)
##
..fastSloc <- function(x, nsamp, k=2, best.r=5, kp, cc, trace=FALSE)
{
    ## NOTES:
    ##  - in the functions rho, psi, and scaledpsi=psi/u (i.e. the weight function)
    ##      is used |x| <= c1
    ##
    ##  - function resdis() to compute the distances is used instead of
    ##      mahalanobis() - slightly faster

    ##  The bisquare rho function:
    ##
    ##              |   x^2/2 - x^4/2*c1^2 + x^6/6*c1^4         |x| <=  c1
    ##  rho(x) =    |
    ##              |   c1^2/6                                  |x| > c1
    ##
    rho <- function(u, cc)
    {
        w <- abs(u) <= cc
        v <- (u^2/2 * (1 - u^2/cc^2 + u^4/(3*cc^4))) * w + (1-w) * (cc^2/6)
        v
    }

    ##  The corresponding psi function: psi = rho'
    ##
    ##              |   x - 2x^3/c1^2 + x^5/c1^4         |x| <=  c1
    ##  psi(x) =    |
    ##              |   0                                |x| > c1
    ##
    ##      using ifelse is 3 times slower
    psi <- function(u, c1)
    {
        ##ifelse(abs(u) < c1, u - 2 * u^3/c1^2 + u^5/c1^4, 0)
        pp <- u - 2 * u^3/c1^2 + u^5/c1^4
        pp*(abs(u) <= c1)
    }

    ## weight function = psi(u)/u
    scaledpsi <- function(u, cc)
    {
        ##ifelse(abs(xx) < c1, xx - 2 * xx^3/c1^2 + xx^5/c1^4, 0)
        pp <- (1 - (u/cc)^2)^2
        pp <- pp * cc^2/6
        pp*(abs(u) <= cc)
    }

    ## the objective function, we solve loss.S(u, s, cc) = b for "s"
    loss.S <- function(u, s, cc)
        mean(rho(u/s, cc))

    norm <- function(x)
        sqrt(sum(x^2))

    ## Returns square root of the mahalanobis distances of x with respect to mu and sigma
    ## Seems to be somewhat more efficient than sqrt(mahalanobis()) - by factor 1.4!
    resdis <- function(x, mu, sigma)
    {
        central <- t(x) - mu
        sqdis <- colSums(solve(sigma, central) * central)
        dis <- sqdis^(0.5)
        dis
    }

    ##  Computes Tukey's biweight objective function (scale)
    ##  (respective to the mahalanobis distances u) using the
    ##  rho() function and the konstants kp and c1
    scaleS <- function(u, kp, c1, initial.sc=median(abs(u))/.6745)
    {
        ## find the scale, full iterations
        maxit <- 200
        eps <- 1e-20
        sc <- initial.sc

        for(i in 1:maxit)
        {
            sc2 <- sqrt(sc^2 * mean(rho(u/sc, c1)) / kp)
            if(abs(sc2/sc - 1) <= eps)
                break
            sc <- sc2
        }
        return(sc)
    }

    ##
    ##  Do "k" iterative reweighting refining steps from "initial.mu, initial.sigma"
    ##
    ##  If "initial.scale" is present, it's used, o/w the MAD is used
    ##
    ##  k = number of refining steps
    ##  conv = 0 means "do k steps and don't check for convergence"
    ##  conv = 1 means "stop when convergence is detected, or the
    ##                 maximum number of iterations is achieved"
    ##  kp and cc = tuning constants of the equation
    ##
    re.s <- function(x, initial.mu, initial.sigma, initial.scale, k, conv, kp, cc)
    {
        n <- nrow(x)
        p <- ncol(x)
        rdis <- resdis(x, initial.mu, initial.sigma)

        if(missing(initial.scale))
        {
            initial.scale <- scale <- median(abs(rdis))/.6745
        } else
        {
            scale <- initial.scale
        }

        ## if conv == 1 then set the max no. of iterations to 50 magic number alert!!!
        if(conv == 1)
            k <- 50

        mu <- initial.mu
        sigma <- initial.sigma
        for(i in 1:k)
        {
            ## do one step of the iterations to solve for the scale
            scale.super.old <- scale
            scale <- sqrt(scale^2 * mean(rho(rdis/scale, cc)) / kp)

            ## now do one step of reweighting with the "improved scale"
            weights <- scaledpsi(rdis/scale, cc)
            W <- weights %*% matrix(rep(1,p), ncol=p)
            xw <- x * W/mean(weights)
            mu.1 <- apply(xw,2,mean)
            res <- x - matrix(rep(1,n),ncol=1) %*% mu.1
            sigma.1 <- t(res) %*% ((weights %*% matrix(rep(1,p), ncol=p)) * res)
            sigma.1 <- (det(sigma.1))^(-1/p) * sigma.1

            if(.isSingular(sigma.1))
            {
                mu.1 <- initial.mu
                sigma.1 <- initial.sigma
                scale <- initial.scale
                break
            }
            if(conv == 1)
            {
                ## check for convergence
                if(norm(mu - mu.1) / norm(mu) < 1e-20)
                    break
                ## magic number alert!!!
            }

            rdis <- resdis(x,mu.1,sigma.1)
            mu <- mu.1
            sigma <- sigma.1
        }

        rdis <- resdis(x,mu,sigma)

        ## get the residuals from the last beta
        return(list(mu.rw = mu.1, sigma.rw=sigma.1, scale.rw = scale))
    }

################################################################################################
    n <- nrow(x)
    p <- ncol(x)

    best.mus <- matrix(0, best.r, p)
    best.sigmas <- matrix(0,best.r*p,p)
    best.scales <- rep(1e20, best.r)
    s.worst <- 1e20
    n.ref <- 1

    for(i in 1:nsamp)
    {
        ## Generate a p+1 subsample in general position and compute
        ##  its mu and sigma. If sigma is singular, generate another
        ##  subsample.
        ##  FIXME ---> the singularity check 'det(sigma) < 1e-7' can
        ##      in some cases give TRUE, e.g. milk although
        ##      .isSingular(sigma) which uses qr(mat)$rank returns
        ##      FALSE. This can result in many retrials and thus
        ##      speeding down!
        ##
        singular <- TRUE
        iii <- 0
        while(singular)
        {
            iii <- iii + 1
            if(iii > 10000)
            {
                cat("\nToo many singular resamples: ", iii, "\n")
                iii <- 0
            }

            indices <- sample(n, p+1)
            xs <- x[indices,]
            mu <- colMeans(xs)
            sigma <- cov(xs)
            ##singular <- det(sigma) < 1e-7
            singular <- .isSingular(sigma)
        }
        sigma <- det(sigma)^(-1/p) * sigma

        ## Perform k steps of iterative reweighting on the elemental set
        if(k > 0)
        {
            ## do the refining
            tmp <- re.s(x=x, initial.mu=mu, initial.sigma=sigma, k=k, conv=0, kp=kp, cc=cc)
            mu.rw <- tmp$mu.rw
            sigma.rw <- tmp$sigma.rw
            scale.rw <- tmp$scale.rw
            rdis.rw <- resdis(x, mu.rw, sigma.rw)
        } else
        {
            ## k = 0 means "no refining"
            mu.rw <- mu
            sigma.rw <- sigma
            rdis.rw <- resdis(x, mu.rw, sigma.rw)
            scale.rw <- median(abs(rdis.rw))/.6745
        }

        if(i > 1)
        {
            ## if this isn't the first iteration....
            ## check whether new mu/sigma belong to the top best results; if so keep
            ## mu and sigma with corresponding scale.
            scale.test <- loss.S(rdis.rw, s.worst, cc)
            if(scale.test < kp)
            {
                s.best <- scaleS(rdis.rw, kp, cc, scale.rw)
                ind <- order(best.scales)[best.r]
                best.scales[ind] <- s.best
                best.mus[ind,] <- mu.rw
                bm1 <- (ind-1)*p;
                best.sigmas[(bm1+1):(bm1+p),] <- sigma.rw
                s.worst <- max(best.scales)
            }
        } else
        {
            ## if this is the first iteration, then this is the best solution anyway...
            best.scales[best.r] <- scaleS(rdis.rw, kp, cc, scale.rw)
            best.mus[best.r,] <- mu.rw
            bm1 <- (best.r-1)*p;
            best.sigmas[(bm1+1):(bm1+p),] <- sigma.rw
        }
    }

    ## do the complete refining step until convergence (conv=1) starting
    ## from the best subsampling candidate (possibly refined)
    super.best.scale <- 1e20
    for(i in best.r:1)
    {
        index <- (i-1)*p
        tmp <- re.s(x=x, initial.mu=best.mus[i,],
                    initial.sigma=best.sigmas[(index+1):(index+p),],
                    initial.scale=best.scales[i],
                    k=0, conv=1, kp=kp, cc=cc)
        if(tmp$scale.rw < super.best.scale)
        {
            super.best.scale <- tmp$scale.rw
            super.best.mu <- tmp$mu.rw
            super.best.sigma <- tmp$sigma.rw
        }
    }

    super.best.sigma <- super.best.scale^2*super.best.sigma
    return(list(
               center=as.vector(super.best.mu),
               cov=super.best.sigma,
               crit=super.best.scale,
               iter=nsamp,
               kp=kp,
               cc=cc,
               method="S-estimates: S-FAST"))
}

##
## Computes S-estimates of multivariate location and scatter by the Ruppert's
##  SURREAL algorithm using Tukey's biweight function
##      data    - the data: a matrix or a data.frame
##      nsamp   - number of random (p+1)-subsamples to be drawn, defaults to 600*p
##      kp, c1  - tunning constants for the  S-estimator with Tukey's biweight
##                function given the breakdown point (bdp) and the dimension p
##
..covSURREAL <- function(data,
                       nsamp,
                       kp,
                       c1,
                       trace=FALSE,
                       tol.inv)
{

## --- now suspended ---
##
## Compute the tunning constant c1 for the  S-estimator with
##  Tukey's biweight function given the breakdown point (bdp)
##  and the dimension p
vcTukey<-function(bdp, p)
{
    cnew <- sqrt(qchisq(1 - bdp, p))
    maxit <- 1000
    epsilon <- 10^(-8)
    error <- 10^6
    contit <- 1
    while(error > epsilon & contit < maxit)
    {
        cold <- cnew
        kp <- vkpTukey(cold,p)
        cnew <- sqrt(6*kp/bdp)
        error <- abs(cold-cnew)
        contit <- contit+1
    }
    if(contit == maxit)
        warning(" Maximum iterations, the last values for c1 and kp are:")
    return (list(c1 = cnew, kp = kp))
}

## Compute the constant kp (used in vcTukey)
vkpTukey<-function(c0, p)
{
   c2 <- c0^2
   intg1 <- (p/2)*pchisq(c2, (p+2))
   intg2 <- (p+2)*p/(2*c2)*pchisq(c2, (p+4))
   intg3 <- ((p+4)*(p+2)*p)/(6*c2^2)*pchisq(c2, (p+6))
   intg4 <- (c2/6)*(1-pchisq(c2, p))
   kp <- intg1-intg2+intg3+intg4
   return(kp)
}

##  Computes Tukey's biweight objective function (scale)
##  (respective to the mahalanobis distances dx)
scaleS <- function(dx, S0, c1, kp, tol)
{
    S <- if(S0 > 0) S0 else median(dx)/qnorm(3/4)
    test <- 2*tol
    rhonew <- NA
    rhoold <- mean(rhobiweight(dx/S,c1)) - kp
    while(test >= tol)
    {
        delta <- rhoold/mean(psibiweight(dx/S,c1)*(dx/(S^2)))
        mmin <- 1
        control <- 0
        while(mmin < 10 & control != 1)
        {
            rhonew <- mean(rhobiweight(dx/(S+delta),c1))-kp
            if(abs(rhonew) < abs(rhoold))
            {
                S <- S+delta
                control <- 1
            }else
            {
                delta <- delta/2
                mmin <- mmin+1
            }
        }
        test <- if(mmin == 10) 0 else (abs(rhoold)-abs(rhonew))/abs(rhonew)
        rhoold <- rhonew
    }
    return(abs(S))
}

rhobiweight <- function(xx, c1)
{
    lessc1 <- (xx^2)/2-(xx^4)/(2*(c1^2))+(xx^6)/(6*c1^4)
    greatc1 <- (c1^2)/6
    lessc1 * abs(xx<c1) + greatc1 * abs(xx>=c1)
}

psibiweight <- function(xx, c1)
{
    (xx - (2*(xx^3))/(c1^2)+(xx^5)/(c1^4))*(abs(xx)<c1)
}

    data <- as.matrix(data)
    n <- nrow(data)
    m <- ncol(data)

    if(missing(nsamp))
        nsamp <- 600*m

    ## TEMPORARY for testing
    ##  ignore kp and c1
    ##v1 <- vcTukey(0.5, m)
    ##if(trace)
    ##    cat("\nSURREAL..: bdp, p, c1, kp=", 0.5, m, v1$c1, v1$kp, "\n")
    ##c1 <- v1$c1
    ##kp <- v1$kp

    ma <- m+1
    mb <- 1/m
    tol <- 1e-5
    Stil <- 10e10

    laux <- 1
    for(l in 1:nsamp)
    {
        ## generate a random subsample of size p+1 and compute its
        ##  mean 'muJ' and covariance matrix 'covJ'
        singular <- TRUE
        while(singular)
        {
            xsetJ <- data[sample(1:n, ma), ]
            muJ <- apply(xsetJ, 2, mean)
            covJ <- var(xsetJ)
            singular <- .isSingular(covJ) || .isSingular(covJ <- det(covJ) ^ -mb * covJ)
        }

        if(l > ceiling(nsamp*0.2))
        {
            if(l == ceiling(nsamp*(0.5)))
                laux <- 2
            if(l == ceiling(nsamp*(0.8)))
                laux <- 4

            epsilon <- runif(1)
            epsilon <- epsilon^laux
            muJ <- epsilon*muJ + (1-epsilon)*mutil
            covJ <- epsilon*covJ + (1-epsilon)*covtil
        }
        covJ <- det(covJ) ^ -mb * covJ

        md<-sqrt(mahalanobis(data, muJ, covJ, tol.inv=tol.inv))
        if(mean(rhobiweight(md/Stil, c1)) < kp)
        {
            if(Stil < 5e10)
                Stil <- scaleS(md, Stil, c1, kp, tol)
            else
                Stil <- scaleS(md, 0, c1, kp, tol)
            mutil <- muJ
            covtil <- covJ
            psi <- psibiweight(md, c1*Stil)
            u <- psi/md
            ubig <- diag(u)
            ux <- ubig%*%data
            muJ <- apply(ux, 2, mean)/mean(u)
            xcenter <- scale(data, center=muJ,scale=FALSE)
            covJ <- t(ubig%*%xcenter)%*%xcenter
            covJ <- det(covJ) ^ -mb * covJ

            control <- 0
            cont <- 1
            while(cont < 3 && control != 1)
            {
                cont <- cont + 1
                md <- sqrt(mahalanobis(data, muJ, covJ, tol.inv=tol.inv))
                if(mean(rhobiweight(md/Stil, c1)) < kp)
                {
                    mutil <- muJ
                    covtil <- covJ
                    control <- 1
                    if(Stil < 5e10)
                        Stil<-scaleS(md, Stil, c1, kp, tol)
                    else
                        Stil<-scaleS(md, 0, c1, kp, tol)
                }else
                {
                    muJ <- (muJ + mutil)/2
                    covJ <- (covJ + covtil)/2
                    covJ <- det(covJ) ^ -mb * covJ
                }
            }
        }
    }

    covtil <- Stil^2 * covtil

    return (list(center = mutil,
                 cov = covtil,
                 crit = Stil,
                 iter = nsamp,
                 kp=kp, cc=c1,
                 method="S-estimates: SURREAL")
            )
}

##
## Compute the S-estimate of location and scatter using the  bisquare function
##
##  NOTE: this function is equivalent to the function ..covSRocke() for computing
##          the rocke type estimates with the following diferences:
##          - the call to .iter.rocke() instead of iter.bic()
##          - the mahalanobis distance returned by iter.bic() is sqrt(mah)
##
##  arp     - not used, i.e. used only in the Rocke type estimates
##  eps     - Iteration convergence criterion
##  maxiter - maximum number of iterations
##  t0, S0  - initial HBDM estimates of center and location
##  nsamp, seed - if not provided T0 and S0, these will be used for CovMve()
##              to compute them
##  initcontrol - if not provided t0 and S0, initcontrol will be used to compute
##              them
##
..covSBic <- function(x, arp, eps, maxiter, t0, S0, nsamp, seed, initcontrol, trace=FALSE){

    dimn <- dimnames(x)
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]
    method <- "S-estimates: bisquare"

    ##  standarization
    mu <- apply(x, 2, median)
    devin <- apply(x, 2, mad)
    x <- scale(x, center = mu, scale=devin)

    ## If not provided initial estimates, compute them as MVE
    ##  Take the raw estimates and standardise the covariance
    ##  matrix to determinant=1
    if(missing(t0) || missing(S0)){
        if(missing(initcontrol))
            init <- CovMve(x, nsamp=nsamp, seed=seed)
        else
            init <- restimate(initcontrol, x)
        cl <- class(init)
        if(cl == "CovMve" || cl == "CovMcd" || cl == "CovOgk"){
            t0 <- init@raw.center
            S0 <- init@raw.cov

##            xinit <- cov.wt(x[init@best,])
##            t0 <- xinit$center
##            S0 <- xinit$cov

            detS0 <-det(S0)
            detS02 <- detS0^(1.0/p)
            S0 <- S0/detS02
        } else{
            t0 <- init@center
            S0 <- init@cov
        }
    }

    ## initial solution
    center <- t0
    cov <- S0

    out <- .iter.bic(x, center, cov, maxiter, eps, trace=trace)
    center <- out$center
    cov <- out$cov
    mah <- out$mah
    iter <- out$iter
    kp <- out$kp
    cc <- out$cc

    mah <- mah^2
    med0 <- median(mah)
    qchis5 <- qchisq(.5,p)
    cov <- (med0/qchis5)*cov
    mah <- (qchis5/med0)*mah

    DEV1 <- diag(devin)
    center  <- DEV1%*%center+mu
    cov <- DEV1%*%cov%*%DEV1
    center <- as.vector(center)

    crit <- determinant(cov, logarithm = FALSE)$modulus[1]
    if(!is.null(nms <- dimn[[2]])) {
        names(center) <- nms
        dimnames(cov) <- list(nms,nms)
    }

    return (list(center=center,
                 cov=cov,
                 crit=crit,
                 method=method,
                 iter=iter,
                 kp=kp,
                 cc=cc,
                 mah=mah)
           )
}

##
## Compute the S-estimate of location and scatter using Rocke type estimator
##
##  see the notes to ..covSBic()
##
..covSRocke <- function(x, arp, eps, maxiter, t0, S0, nsamp, seed, initcontrol, trace=FALSE){

    dimn <- dimnames(x)
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]
    method <- "S-estimates: Rocke type"

    ##  standarization
    mu <- apply(x, 2, median)
    devin <- apply(x, 2, mad)
    x <- scale(x, center = mu, scale=devin)

    ## If not provided initial estimates, compute them as MVE
    ##  Take the raw estimates and standardise the covariance
    ##  matrix to determinant=1
    if(missing(t0) || missing(S0)){
        if(missing(initcontrol))
            init <- CovMve(x, nsamp=nsamp, seed=seed)
        else
            init <- restimate(initcontrol, x)
        cl <- class(init)
        if(cl == "CovMve" || cl == "CovMcd" || cl == "CovOgk"){
            t0 <- init@raw.center
            S0 <- init@raw.cov

##            xinit <- cov.wt(x[init@best,])
##            t0 <- xinit$center
##            S0 <- xinit$cov

            detS0 <-det(S0)
            detS02 <- detS0^(1.0/p)
            S0 <- S0/detS02
        } else{
            t0 <- init@center
            S0 <- init@cov
        }
    }

    ## initial solution
    center <- t0
    cov <- S0

    out <- .iter.rocke(x, center, cov, maxiter, eps, alpha=arp, trace=trace)
    center <- out$center
    cov <- out$cov
    mah <- out$mah
    iter <- out$iter

#    mah <- mah^2
    med0 <- median(mah)
    qchis5 <- qchisq(.5,p)
    cov <- (med0/qchis5)*cov
    mah <- (qchis5/med0)*mah

    DEV1 <- diag(devin)
    center  <- DEV1%*%center+mu
    cov <- DEV1%*%cov%*%DEV1
    center <- as.vector(center)

    crit <- determinant(cov, logarithm = FALSE)$modulus[1]
    if(!is.null(nms <- dimn[[2]])) {
        names(center) <- nms
        dimnames(cov) <- list(nms,nms)
    }

    return (list(center=center,
                 cov=cov,
                 crit=crit,
                 method=method,
                 iter=iter,
                 kp=0,
                 cc=0,
                 mah=mah)
           )
}

.iter.bic <- function(x, center, cov, maxiter, eps, trace)
{

    stepS <- function(x, b, cc, mah, n, p)
    {
        p2 <- (-1/p)
        w <- w.bi(mah,cc)
        w.v <- matrix(w,n,p)
        tn <- colSums(w.v*x)/sum(w)
        xc <- x-t(matrix(tn,p,n))
        vn <- t(xc)%*%(xc*w.v)
        eii <- eigen(vn, symmetric=TRUE, only.values = TRUE)
        eii <- eii$values
        aa <- min(eii)
        bb <- max(eii)
        condn <- bb/aa
        cn <- NULL
        fg <- TRUE

        if(condn > 1.e-12)
        {
            fg <- FALSE
            cn <- (det(vn)^p2)*vn
        }

        out1 <- list(tn, cn, fg)
        names(out1) <- c("center", "cov", "flag")
        return(out1)
    }

    scale.full <- function(u, b, cc)
    {
        #find the scale - full iterations
        max.it <- 200
        s0 <- median(abs(u))/.6745
        s.old <- s0
        it <- 0
        eps <- 1e-4
        err <- eps + 1
        while(it < max.it && err > eps)
        {
            it <- it+1
            s.new <- s.iter(u,b,cc,s.old)
            err <- abs(s.new/s.old-1)
            s.old <- s.new
        }
        return(s.old)
    }

    rho.bi <- function(x, cc)
    {
        ## bisquared rho function
        dif <- abs(x) - cc < 0
        ro1 <- (x^2)/2 - (x^4)/(2*(cc^2)) + (x^6)/(6*(cc^4))
        ro1[!dif] <- cc^2/6
        return(ro1)
    }

    psi.bi <- function(x, cc)
    {
        ## bisquared psi function
        dif <- abs(x) - cc < 0
        psi1 <- x - (2*(x^3))/(cc^2) + (x^5)/(cc^4)
        psi1[!dif] <- 0
        return(psi1)
    }

    w.bi <- function(x, cc)
    {
        ## bicuadratic w function
        dif <- abs(x) - cc < 0
        w1 <- 1 - (2*(x^2))/(cc^2) + (x^4)/(cc^4)
        w1[!dif] <- 0
        w1
    }

    s.iter <- function(u, b, cc, s0)
    {
        # does one step for the scale
        ds <- u/s0
        ro <- rho.bi(ds, cc)
        snew <- sqrt(((s0^2)/b) * mean(ro))
        return(snew)
    }

##################################################################################

    ##  main routine for .iter.bic()
    cc <- 1.56
    b <- cc^2/12

    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]

    if(trace)
        cat("\nbisquare.: bdp, p, c1, kp=", 0.5, p, cc, b, "\n")

    mah <- sqrt(mahalanobis(x, center, cov))
    s0 <- scale.full(mah, b, cc)
    mahst <- mah/s0

    ## main loop
    u <- eps + 1
    fg <- FALSE
    it <- 0
    while(u > eps & !fg & it < maxiter)
    {
        it <- it + 1
        out <- stepS(x, b, cc, mahst, n, p)
        fg <- out$flag
        if(!fg)
        {
            center <- out$center
            cov <- out$cov
            mah <- sqrt(mahalanobis(x, center, cov))
            s <- scale.full(mah, b, cc)
            mahst <- mah/s
            u <- abs(s-s0)/s0
            s0 <- s
        }
    }

    list(center=center, cov=cov, mah=mah, iter=it, kp=b, cc=cc)
}


.iter.rocke <- function(x, mu, v, maxit, tol, alpha, trace)
{


rho.rk2 <- function(t, p, alpha)
{
    x <- t
    z  <- qchisq(1-alpha, p)
    g <- min(z/p - 1, 1)
    uu1 <- (x <= 1-g)
    uu2 <- (x > 1-g & x <= 1+g)
    uu3 <- (x > 1+g)
    zz <- x
    x1 <- x[uu2]
    dd <- ((x1-1)/(4*g))*(3-((x1-1)/g)^2)+.5
    zz[uu1] <- 0
    zz[uu2] <- dd
    zz[uu3] <- 1
    return(zz)
}

w.rk2 <- function(t, p, alpha)
{
    x <- t
    z <- qchisq(1-alpha, p)
    g <- min(z/p - 1, 1)
    uu1 <- (x <= (1-g) )
    uu2 <- (x > 1-g & x <= 1+g)
    uu3 <- (x > (1+g))
    zz <- x
    x1 <- x[uu2]
    dd <- (3/(4*g))*(1-((x1-1)/g)^2)
    zz[uu1] <- 0
    zz[uu2] <- dd
    zz[uu3] <- 0
    return(zz)
}

disesca <- function(x, mu, v, alpha, delta)
{
    ##  calculate square roots of Mahalanobis distances and scale
    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]
    p2 <- (1/p)
    p3 <- -p2
    maha <- mahalanobisfl(x,mu,v)
    dis <- maha$dist
    fg <- maha$flag
    detv <- maha$det

    sig <- NULL
    v1 <- NULL
    dis1 <- NULL
    if(!fg)
    {
        dis1 <- dis*(detv^p2)
        v1 <- v*(detv^p3)
        sig <- escala(p, dis1, alpha, delta)$sig
    }
    out <- list(dis=dis1, sig=sig, flag=fg, v=v1)
    return(out)
}

fespro <- function(sig, fp, fdis, falpha, fdelta)
{
    z <- fdis/sig
    sfun <- mean(rho.rk2(z, fp, falpha)) - fdelta
    return(sfun)
}

escala <- function(p, dis, alpha, delta)
{
    ##  find initial interval for scale
    sig <- median(dis)
    a <- sig/100
    ff <- fespro(a, p, dis, alpha, delta)
    while(ff < 0) 
    {
        a <- a/10
        ff <- fespro(a, p, dis, alpha, delta)
    }

    b <- sig*10
    ff <- fespro(b, p, dis, alpha, delta)
    while(ff>0) 
    {
        b <- b*10
        ff <- fespro(b, p, dis, alpha, delta)
    }

    ##  find scale (root of fespro)
    sig.res <- uniroot(fespro, lower=a, upper=b, fp=p, fdis=dis, falpha=alpha, fdelta=delta)
    sig <- sig.res$root
    mens <- sig.res$message
    ans <- list(sig=sig, message=mens)
    return(ans)
}

gradien <- function(x, mu1, v1, mu0, v0, sig0, dis0, sig1, dis1, alpha)
{
    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]
    delta <- 0.5*(1-p/n)
    decre <- 0.
    mumin <- mu0
    vmin <- v0
    dismin <- dis0
    sigmin <- sig0
    for(i in 1:10) {
        gamma <- i/10
        mu <- (1-gamma)*mu0+gamma*mu1
        v <- (1-gamma)*v0+gamma*v1
        dis.res <- disesca(x,mu,v,alpha,delta)
        if(!dis.res$flag) {
            sigma <- dis.res$sig
            decre <- (sig0-sigma)/sig0
            if(sigma < sigmin) {
                mumin <- mu
                vmin <- dis.res$v
                sigmin <- dis.res$sig
                dismin <- dis.res$dis
            }
        }
    }
    ans <- list(mu=mumin, v=vmin, sig=sigmin, dis=dismin)
    return(ans)
}



descen1 <- function(x, mu, v, alpha, sig0, dis0)
{
    ##  one step of reweighted
    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]
    p2 <- -(1/p)
    delta <- 0.5*(1-p/n)
    daux  <- dis0/sig0
    w <- w.rk2(daux, p, alpha)
    w.v <- w %*% matrix(1,1,p)
    mui <- colSums(w.v * x)/sum(w)
    xc <- x- (matrix(1,n,1) %*% mui)
    va <- t(xc) %*% (xc*w.v)
    disi.res <- disesca(x,mui,va,alpha,delta)
    if(disi.res$flag)
    {
        disi.res$sig <- sig0+1
        disi.res$v <- va
    }
    disi <- disi.res$dis
    sigi <- disi.res$sig
    flagi <- disi.res$flag
    deti<-disi.res$det
    vi <- disi.res$v

    ans <- list(mui=mui, vi=vi, sig=sigi, dis=disi, flag=flagi)
    return(ans)
}

invfg <- function(x, tol=1.e-12)
{
    dm <- dim(x)
    p <- dm[1]
    y <- svd(x)
    u <- y$u
    v <- y$v
    d <- y$d
    cond <- min(d)/max(d)
    fl <- FALSE
    if(cond < tol)
        fl <- TRUE

    det <- prod(d)
    invx <- NULL
    if(!fl)
        invx <- v%*%diag(1/d)%*%t(u)
    out <- list(inv=invx, flag=fl, det=det)
    return(out)
}


mahalanobisfl <- function(x, center, cov, tol=1.e-12)
{
    invcov <- invfg(cov,tol)
    covinv <- invcov$inv
    fg <- invcov$flag
    det1 <- invcov$det
    dist <- NULL
    if(!fg)
        dist <- mahalanobis(x, center, covinv, inverted=TRUE)

    out <- list(dist=dist, flag=fg, det=det1)
    return(out)
}

##################################################################################

    ##  main routine for .iter.rocke()
    if(trace)
    {
        cat("\nWhich initial t0 and S0 are we using: \n")
        print(mu)
        print(v)
    }
    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]
    delta <- 0.5*(1-p/n)
    dis.res <- disesca (x, mu, v, alpha, delta)
    vmin <- v
    mumin <- mu
    sigmin <- dis.res$sig
    dismin <- dis.res$dis
    it <- 1
    decre <- tol+1
    while(decre > tol && it <= maxit)
    {
        ##  do reweighting
        des1 <- descen1(x, mumin, vmin, alpha, sigmin, dismin)
        mu1 <- des1$mui
        v1 <- des1$vi
        dis1 <- des1$dis
        sig1 <- des1$sig
        decre1 <- (sigmin-sig1)/sigmin
        if(trace)
        {
            cat("\nit, decre1:",it,decre1, "\n")
            print(des1)
        }
        if(decre1 <= tol) 
        {
            grad <- gradien(x, mu1, v1, mumin, vmin, sigmin, dismin, sig1, dis1, alpha)
            mu2 <- grad$mu
            v2 <- grad$v
            sig2 <- grad$sig
            dis2<-grad$dis
            decre2 <- (sigmin-sig2)/sigmin
            decre <- max(decre1,decre2)
            if(decre1 >= decre2) {
                mumin <- mu1
                vmin <- v1
                dismin <- dis1
                sigmin <- sig1
            }
            if(decre1 < decre2) {
                mumin <- mu2
                vmin <- v2
                dismin <- dis2
                sigmin <- sig2
            }
        }
        
        if(decre1 > tol) 
        {
            decre <- decre1
            mumin <- mu1
            vmin <- v1
            dismin <- dis1
            sigmin <- sig1
        }
        it <- it+1
    }
    ans <- list(center=mumin, cov=vmin, scale=sigmin, mah=dismin, iter=it)
    return(ans)
}
