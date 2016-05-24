
norm <- function(x) sqrt(sum(x^2))

########## REFERENCES ##########
#####
##### Fletcher, R. (1987)
##### Practical Methods of Optimization, second edition.
##### John Wiley, Chichester.
#####
##### Nocedal, J. and Wright, S. J. (1999)
##### Numerical Optimization.
##### Springer-Verlag, New York.
#####
##### See Section 5.1 of Fletcher
##### See Section 4.2 of Nocedal and Wright
#####
################################

########## COMMENT ##########
##### Our method using one eigendecomposition per iteration is not fastest.
##### Both books recommend using multiple Cholesky decompositions instead.
##### But the eigendecomposition method is simpler to program, easier to
##### understand (which is why both books use it for their theoretical
##### explanation), and hopefully more bulletproof.
#####
##### Our idea for this comes from the way mvrnorm in the MASS package also
##### uses eigendecomposition rather than Cholesky -- also because bulletproof
##### is better than fast.
#############################

trust <- function(objfun, parinit, rinit, rmax, parscale,
    iterlim = 100, fterm = sqrt(.Machine$double.eps),
    mterm = sqrt(.Machine$double.eps),
    minimize = TRUE, blather = FALSE, ...)
{
    if (! is.numeric(parinit))
       stop("parinit not numeric")
    if (! all(is.finite(parinit)))
       stop("parinit not all finite")
    d <- length(parinit)
    if (missing(parscale)) {
        rescale <- FALSE
    } else {
        rescale <- TRUE
        if (length(parscale) != d)
           stop("parscale and parinit not same length")
        if (! all(parscale > 0))
           stop("parscale not all positive")
        if (! all(is.finite(parscale) & is.finite(1 / parscale)))
           stop("parscale or 1 / parscale not all finite")
    }
    if (! is.logical(minimize))
       stop("minimize not logical")

    r <- rinit
    theta <- parinit
    out <- try(objfun(theta, ...))
    if (inherits(out, "try-error")) {
        warning("error in first call to objfun")
        return(list(error = out, argument = theta, converged = FALSE,
            iterations = 0))
    }
    check.objfun.output(out, minimize, d)
    if (! is.finite(out$value))
        stop("parinit not feasible")
    accept <- TRUE

    if (blather) {
        theta.blather <- NULL
        theta.try.blather <- NULL
        type.blather <- NULL
        accept.blather <- NULL
        r.blather <- NULL
        stepnorm.blather <- NULL
        rho.blather <- NULL
        val.blather <- NULL
        val.try.blather <- NULL
        preddiff.blather <- NULL
    }

    for (iiter in 1:iterlim) {

        if (blather) {
            theta.blather <- rbind(theta.blather, theta)
            r.blather <- c(r.blather, r)
            if (accept)
                val.blather <- c(val.blather, out$value)
            else
                val.blather <- c(val.blather, out.value.save)
        }

        if (accept) {
            B <- out$hessian
            g <- out$gradient
            f <- out$value
            out.value.save <- f
            if (rescale) { 
                B <- B / outer(parscale, parscale)
                g <- g / parscale
            }
            if (! minimize) {
                B <- (- B)
                g <- (- g)
                f <- (- f)
            }
            eout <- eigen(B, symmetric = TRUE)
            gq <- as.numeric(t(eout$vectors) %*% g)
        }

        ########## solve trust region subproblem ##########

        ##### try for Newton #####
        is.newton <- FALSE
        if (all(eout$values > 0)) {
            ptry <- as.numeric(- eout$vectors %*% (gq / eout$values))
            if (norm(ptry) <= r)
                is.newton <- TRUE
        }

        ##### non-Newton #####
        if (! is.newton) {
            lambda.min <- min(eout$values)
            beta <- eout$values - lambda.min
            imin <- beta == 0
            C1 <- sum((gq / beta)[! imin]^2)
            C2 <- sum(gq[imin]^2)
            C3 <- sum(gq^2)
            if (C2 > 0 || C1 > r^2) {
                is.easy <- TRUE
                is.hard <- (C2 == 0)
                ##### easy cases #####
                beta.dn <- sqrt(C2) / r
                beta.up <- sqrt(C3) / r
                fred <- function(beep) {
                    if (beep == 0) {
                        if (C2 > 0)
                            return(- 1 / r)
                        else
                            return(sqrt(1 / C1) - 1 / r)
                    }
                    return(sqrt(1 / sum((gq / (beta + beep))^2)) - 1 / r)
                }
                if (fred(beta.up) <= 0) {
                    uout <- list(root = beta.up)
                } else if (fred(beta.dn) >= 0) {
                    uout <- list(root = beta.dn)
                } else {
                    uout <- uniroot(fred, c(beta.dn, beta.up))
                }
                wtry <- gq / (beta + uout$root)
                ptry <- as.numeric(- eout$vectors %*% wtry)
            } else {
                is.hard <- TRUE
                is.easy <- FALSE
                ##### hard-hard case #####
                wtry <- gq / beta
                wtry[imin] <- 0
                ptry <- as.numeric(- eout$vectors %*% wtry)
                utry <- sqrt(r^2 - sum(ptry^2))
                if (utry > 0) {
                    vtry <- eout$vectors[ , imin, drop = FALSE]
                    vtry <- vtry[ , 1]
                    ptry <- ptry + utry * vtry
                }
            }
        }

        ########## predicted versus actual change ##########
        preddiff <- sum(ptry * (g + as.numeric(B %*% ptry) / 2))
        if (rescale) {
            theta.try <- theta + ptry / parscale
        } else {
            theta.try <- theta + ptry
        }
        out <- try(objfun(theta.try, ...))
        if (inherits(out, "try-error"))
            break
        check.objfun.output(out, minimize, d)
        ftry <- out$value
        if (! minimize)
            ftry <- (- ftry)
        rho <- (ftry - f) / preddiff

        ########## termination test ##########
        if (ftry < Inf) {
            is.terminate <- abs(ftry - f) < fterm || abs(preddiff) < mterm
        } else {
            is.terminate <- FALSE
            rho <- (- Inf)
        }

        ##### adjustments #####
        if (is.terminate) {
            if (ftry < f) {
                accept <- TRUE
                theta <- theta.try
            }
        } else {
            if (rho < 1 / 4) {
                accept <- FALSE
                r <- r / 4
            } else {
                accept <- TRUE
                theta <- theta.try
                if (rho > 3 / 4 && (! is.newton))
                    r <- min(2 * r, rmax)
            }
        }

        if (blather) {
            theta.try.blather <- rbind(theta.try.blather, theta.try)
            val.try.blather <- c(val.try.blather, out$value)
            accept.blather <- c(accept.blather, accept)
            preddiff.blather <- c(preddiff.blather, preddiff)
            stepnorm.blather <- c(stepnorm.blather, norm(ptry))
            if (is.newton) {
                mytype <- "Newton"
            } else {
                if (is.hard) {
                    if (is.easy) {
                        mytype <- "hard-easy"
                    } else {
                        mytype <- "hard-hard"
                    }
                } else {
                    mytype <- "easy-easy"
                }
            }
            type.blather <- c(type.blather, mytype)
            rho.blather <- c(rho.blather, rho)
        }

        if (is.terminate)
            break
    }

    if (inherits(out, "try-error")) {
        out <- list(error = out, argument = theta.try, converged = FALSE)
    } else {
        out <- try(objfun(theta, ...))
        if (inherits(out, "try-error")) {
            out <- list(error = out)
            warning("error in last call to objfun")
        } else {
            check.objfun.output(out, minimize, d)
        }
        out$argument <- theta
        out$converged <- is.terminate
    }
    out$iterations <- iiter
    if (blather) {
        dimnames(theta.blather) <- NULL
        out$argpath <- theta.blather
        dimnames(theta.try.blather) <- NULL
        out$argtry <- theta.try.blather
        out$steptype <- type.blather
        out$accept <- accept.blather
        out$r <- r.blather
        out$rho <- rho.blather
        out$valpath <- val.blather
        out$valtry <- val.try.blather
        if (! minimize)
            preddiff.blather <- (- preddiff.blather)
        out$preddiff <- preddiff.blather
        out$stepnorm <- stepnorm.blather
    }
    return(out)
}

check.objfun.output <- function(obj, minimize, dimen)
{
    if (! is.list(obj))
        stop("objfun returned object that is not a list")
    foo <- obj$value
    if (is.null(foo))
        stop("objfun returned list that does not have a component 'value'")
    if (! is.numeric(foo))
        stop("objfun returned value that is not numeric")
    if (length(foo) != 1)
        stop("objfun returned value that is not scalar")
    if (is.na(foo) || is.nan(foo))
        stop("objfun returned value that is NA or NaN")
    if (minimize && foo == (-Inf))
        stop("objfun returned -Inf value in minimization")
    if ((! minimize) && foo == Inf)
        stop("objfun returned +Inf value in maximization")
    if (is.finite(foo)) {
        bar <- obj$gradient
        if (is.null(bar))
            stop("objfun returned list without component 'gradient' when value is finite")
        if (! is.numeric(bar))
            stop("objfun returned gradient that is not numeric")
        if (length(bar) != dimen)
            stop(paste("objfun returned gradient that is not vector of length", dimen))
        if (! all(is.finite(bar)))
            stop("objfun returned gradient not having all elements finite")
        baz <- obj$hessian
        if (is.null(baz))
            stop("objfun returned list without component 'hessian' when value is finite")
        if (! is.numeric(baz))
            stop("objfun returned hessian that is not numeric")
        if (! is.matrix(baz))
            stop("objfun returned hessian that is not matrix")
        if (! all(dim(baz) == dimen))
            stop(paste("objfun returned hessian that is not", dimen, "by", dimen, "matrix"))
        if (! all(is.finite(baz)))
            stop("objfun returned hessian not having all elements finite")
    }
    return(TRUE)
}
