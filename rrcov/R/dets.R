##
## Input:
##      x      - a data matrix of size (n,p)
##      hsets.init
##      scalefn
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

..detSloc <- function(x,
                      hsets.init=NULL,
                      save.hsets=missing(hsets.init), full.h=save.hsets,
                      scalefn,
                      maxisteps = 200, warn.nonconv.csteps = FALSE,
                      ## k=0,              # no need of preliminary refinement
                      ## best.r=6,         # iterate till convergence on all 6 sets
                      kp,
                      cc,
                      trace=as.integer(trace))

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
    re.s <- function(x, initial.mu, initial.sigma, initial.scale, k, conv, maxisteps=200, kp, cc)
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
            k <- maxisteps

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

            if(i >= k & warn.nonconv.csteps)
                warning("Convergence not reached up to max.number of iterations: ", k)
        }

        iter <- i
        rdis <- resdis(x,mu,sigma)

        ## get the residuals from the last beta
        ## return the number of iterations
        return(list(mu.rw = mu.1, sigma.rw=sigma.1, scale.rw = scale, iter=iter))
    }

################################################################################################
    n <- nrow(x)
    p <- ncol(x)
    h <- h.alpha.n(0.5, n, p)

    ## Center and scale the data
    z <- doScale(x, center=median, scale=scalefn)
    z.center <- z$center
    z.scale <- z$scale
    z <- z$x

    ## Assume that 'hsets.init' already contains h-subsets: the first h observations each
    if(is.null(hsets.init))
    {
    	hsets.init <- r6pack(z, h=h, full.h=full.h, scaled=TRUE, scalefn=scalefn)
    	dh <- dim(hsets.init)
    } else { ## user specified, (even just *one* vector):
    	if(is.vector(hsets.init)) hsets.init <- as.matrix(hsets.init)
    	dh <- dim(hsets.init)
    	if(!is.matrix(hsets.init) || dh[1] < h || dh[2] < 1)
    	    stop("'hsets.init' must be a  h' x L  matrix (h' >= h) of observation indices")
    	if(full.h && dh[1] != n)
    	    stop("When 'full.h' is true, user specified 'hsets.init' must have n rows")
    }

    ## Construction of h-subset: take the first h observations
    hsets <- hsets.init[1:h, , drop=FALSE]       ## select the first h ranked observations
    nsets <- NCOL(hsets)
    hset.csteps <- integer(nsets)

    best.r <- nsets
    k <- 0
    best.mus <- matrix(0, best.r, p)
    best.sigmas <- matrix(0,best.r*p,p)
    best.scales <- rep(1e20, best.r)
    s.worst <- 1e20
    n.ref <- 1

    ## Iterate now through the initial h-subsets
    for(i in 1:nsets)
    {

        if(trace) {
            if(trace >= 2)
                cat(sprintf("H-subset %d = observations c(%s):\n-----------\n",
                            i, paste(hsets.init[1:h,i], collapse=", ")))
            else
                cat(sprintf("H-subset %d: ", i))
        }

        indices <- hsets[, i]     # start with the i-th initial set
        xs <- x[indices,]
        mu <- colMeans(xs)
        sigma <- cov(xs)
        singular <- .isSingular(sigma)

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
                    k=0, conv=1, maxisteps=maxisteps, kp=kp, cc=cc)
        if(tmp$scale.rw < super.best.scale)
        {
            super.best.scale <- tmp$scale.rw
            super.best.mu <- tmp$mu.rw
            super.best.sigma <- tmp$sigma.rw
            ind.best <- i # to determine which subset gives best results.
        }

        hset.csteps[i] <- tmp$iter # how many I-steps necessary to converge.
        if(trace) cat(sprintf("%3d csteps, scale=%g", tmp$iter, tmp$scale.rw))
    }

    super.best.sigma <- super.best.scale^2*super.best.sigma
    return(list(
               center=as.vector(super.best.mu),
               cov=super.best.sigma,
               crit=super.best.scale,
               iter=hset.csteps[ind.best],
               iBest = ind.best,
               hset.csteps = hset.csteps,
               initHsets=if(save.hsets) hsets.init,
               kp=kp,
               cc=cc,
               method="S-estimates: DET-S"))
}
