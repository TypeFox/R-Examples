#######################################################
## Main methods                                      ##
#######################################################

####
## laplace approximated log likelihood
## use the following "algorithm":
##
## 1. find beta, b for given values of theta
## 2. find theta for given values of beta, b
## 3. iterate
####

##' Robust estimation of linear mixed effects models, for hierarchical
##' nested and non-nested, e.g., crossed, datasets.
##'
##' \describe{
##' \item{Overview:}{
##' This function implements a robust approach of fitting linear mixed
##' effect models. It can be used much like the function
##' \code{\link[lme4]{lmer}} in the package \code{lme4}. The supported
##' models are the same as for \code{\link[lme4]{lmer}} (gaussian
##' family only). The robust approach used is based on the
##' robustification of the scoring equations and an application of the
##' Design Adaptive Scale approach.
##'
##' Example analyses and theoretical details on the method
##' are available in the vignette (see \code{vignette("rlmer")}).
##'
##' Models are specified using the \code{formula} argument, using the
##' same syntax as for \code{\link[lme4]{lmer}}. Additionally, one
##' also needs to specify what robust scoring or weight functions are
##' to be used (arguments starting with \code{rho.}). By default a
##' smoothed version of the Huber function is used. Furthermore, the
##' \code{method} argument can be used to speed up computations at the
##' expense of accuracy of the results.
##' }
##' 
##' \item{Computation methods:}{
##' Currently, there are two different methods available for fitting
##' models. They only differ in how the consistency factors for the
##' Design Adaptive Scale estimates are computed.
##' Available fitting methods for theta and sigma.e:
##' \itemize{
##' \item \code{DAStau} (default): For this method, the consistency
##'         factors are computed using numerical quadrature. This is
##'         slower but yields more accurate results. This is the direct
##'         analogue to the DAS-estimate in robust linear regression.
##' \item \code{DASvar}:
##'         This method computes the consistency factors using a
##'         direct approximation which is faster but less accurate.
##'         For complex models with correlated random effects with
##'         more than one correlation term, this is the only
##'         method available.
##' }
##' }
##' 
##' \item{Weight functions:}{
##' The tuning parameters of the weight functions \dQuote{rho} can be
##' used to adjust robustness and efficiency of the resulting
##' estimates (arguments \code{rho.e}, \code{rho.b},
##' \code{rho.sigma.e} and \code{rho.sigma.b}). Better robustness will
##' lead to a decrease of the efficiency. By default, the tuning
##' parameters are set to yield estimates with approximately 95\%
##' efficiency for the fixed effects. The variance components are
##' estimated with a lower efficiency but better robustness properties.
##'
##' One has to use different weight functions and tuning parameters
##' for simple variance components and for such including correlation
##' parameters. By default, they are chosen appropriately to the model
##' at hand. However, when using the \code{rho.sigma.e} and
##' \code{rho.sigma.b} arguments, it is up to the used to specify
##' the appropriate function.
##' \itemize{
##' \item For simple variance components and the residual error scale
##' use the function \code{\link{psi2propII}} to change the tuning
##' parameters. The is similar to Proposal II in the location-scale
##' problem (i.e., using the squared robustness weights of the
##' location estimate for the scale estimate; otherwise the scale
##' estimate is not robust).
##'
##' \item For random effects modeled with correlation parameters
##' (referred to as nondiagonal case below), use the
##' \code{\link{chgDefaults}} function to change the tuning
##' parameters. The parameter estimation problem is multivariate,
##' unlike the case without correlation where the problem was
##' univariate. For the employed estimator, this amounts to switching
##' from simple scale estimates to estimating correlation
##' matrices. Therefore different weight functions have to be
##' used. Squaring of the weights (using the function
##' \code{\link{psi2propII}}) is no longer necessary. To yield
##' estimates with the same efficiency, the tuning parameters for the
##' nondiagonal are generally larger than for the simple case. As a
##' rule of thumb, one may use the squared tuning parameters of the
##' simple case for the nondiagonal case.  }
##'
##' Tables of tuning factors are given in the vignette
##' (\code{vignette("rlmer")}). For the smoothed Huber function the
##' tuning parameters to get approximately 95\% efficiency are
##' \eqn{k=2.28}{k=2.28} for simple variance components and
##' \eqn{k=5.11}{k=5.11} for variance components including correlation
##' parameters.
##' }
##'
##' \item{Specifying (multiple) weight functions:}{
##' If custom weight functions are specified using the argument
##' \code{rho.b} (\code{rho.e}) but the argument \code{rho.sigma.b}
##' (\code{rho.sigma.e}) is missing, then the squared weights are used
##' for simple variance components and the regular weights are used for
##' variance components including correlation parameters. The same
##' tuning parameters will be used, to get higher efficiency one has
##' to specify the tuning parameters by hand using the
##' \code{\link{psi2propII}} and \code{\link{chgDefaults}} functions.
##'
##' To specify separate weight functions \code{rho.b} and
##' \code{rho.sigma.b} for different variance components, it is
##' possible to pass a list instead of a psi_func object. The list
##' entries correspond to the groups as shown by \code{VarCorr(.)}
##' when applied to the model fitted with \code{lmer}. A set of
##' correlated random effects count as just one group.
##' }
##' }
##'
##' @title Robust linear mixed models
##' @param formula a two-sided linear formula object describing the
##'   fixed-effects part of the model, with the response on the left of
##'   a \code{~} operator and the terms, separated by \code{+}
##'   operators, on the right.  The vertical bar character \code{"|"}
##'   separates an expression for a model matrix and a grouping factor.
##' @param data an optional data frame containing the variables named
##'   in \code{formula}.  By default the variables are taken from the
##'   environment from which \code{lmer} is called.
##' @param ... Additional parameters passed to lmer to find the
##'   initial estimates. See \code{\link[lme4]{lmer}}.
##' @param method method to be used for estimation of theta and sigma,
##'   see Details.
##' @param rho.e object of class psi_func, specifying the functions to
##'   use for the huberization of the residuals.
##' @param rho.b object of class psi_func or list of such objects
##'   (see Details), specifying the functions to use for the
##'   huberization of the random effects.
##' @param rho.sigma.e object of class psi_func, specifying the
##'   weight functions to use for the huberization of the residuals when
##'   estimating the variance components, use the
##'   \code{\link{psi2propII}} function to specify squared weights
##'   and custom tuning parameters.
##' @param rho.sigma.b (optional) object of class psi_func or list of
##'   such objects, specifying the weight functions to use for the
##'   huberization of the random effects when estimating the variance
##'   components (see Details). Use \code{\link{psi2propII}} to specify
##'   squared weights and custom tuning parameters or
##'   \code{\link{chgDefaults}} for regular weights for variance components
##'   including correlation parameters.
##' @param rel.tol relative tolerance used as criteria in the fitting
##'   process.
##' @param max.iter maximum number of iterations allowed.
##' @param verbose verbosity of output. Ranges from 0 (none) to 3
##'   (a lot of output)
##' @param doFit logical scalar. When \code{doFit = FALSE} the model
##'   is not fit but instead a structure with the model matrices for the
##'   random-effects terms is returned (used to speed up tests). When
##'   \code{doFit = TRUE}, the default, the model is fit immediately.
##' @param init optional lmerMod- or rlmerMod-object to use for starting
##'   values, a list with elements \sQuote{fixef}, \sQuote{u},
##'   \sQuote{sigma}, \sQuote{theta}, or a function producing an lmerMod
##'   object.
##' @return object of class rlmerMod.
##' @seealso \code{\link[lme4]{lmer}}
##' @author Manuel Koller, with thanks to Vanda Lourenço for improvements.
##' @keywords models
##' @examples
##' ## dropping of VC
##' system.time(print(rlmer(Yield ~ (1|Batch), Dyestuff2, method="DASvar")))
##'
##' \dontrun{
##'   ## Default method "DAStau"
##'   system.time(rfm.DAStau <- rlmer(Yield ~ (1|Batch), Dyestuff))
##'   summary(rfm.DAStau)
##'   ## DASvar method (faster, less accurate)
##'   system.time(rfm.DASvar <- rlmer(Yield ~ (1|Batch), Dyestuff,
##'                                   method="DASvar"))
##'   ## compare the two
##'   compare(rfm.DAStau, rfm.DASvar)
##' 
##'   ## Fit variance components with higher efficiency
##'   ## psi2propII yields squared weights to get robust estimates
##'   rlmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin,
##'         rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
##'         rho.sigma.b = psi2propII(smoothPsi, k = 2.28))
##'
##'   ## use chgDefaults for variance components including
##'   ## correlation terms (regular, non squared weights suffice)
##'   rlmer(Reaction ~ Days + (Days|Subject), sleepstudy,
##'         rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
##'         rho.sigma.b = chgDefaults(smoothPsi, k = 5.11, s=10))
## }
##'
##' @importFrom lme4 lmer
##' @importFrom stats getCall
##' @export
rlmer <- function(formula, data, ..., method = "DAStau",
                  rho.e = smoothPsi, rho.b = smoothPsi,
                  rho.sigma.e, rho.sigma.b, rel.tol = 1e-8,
                  max.iter = 40*(r+1)^2, verbose = 0,
                  doFit = TRUE, init)
{
    lcall <- match.call()
    if (missing(init) || is.null(init) || is.list(init)) {
        lcall2 <- lcall
        lcall2[setdiff(names(formals(rlmer)), names(formals(lmer)))] <- NULL
        lcall2$doFit <- NULL
        lcall2$REML <- TRUE
        lcall2[[1]] <- as.name("lmer")
        pf <- parent.frame()
        linit <- eval(lcall2, pf)
        if (!missing(init) && is.list(init)) {
            ## check sanity of input
            stopifnot(length(init$fixef) == length(fixef(linit)),
                      length(init$u) == length(getME(linit, "u")),
                      length(init$sigma) == length(sigma(linit)),
                      length(init$theta) == length(getME(linit, "theta")))
            ## convert object to rlmerMod
            linit <- as(linit, "rlmerMod")
            ## set all of the initial parameters, but do not fit yet
            setFixef(linit, unname(init$fixef))
            setSigma(linit, init$sigma)
            setTheta(linit, init$theta, fit.effects=FALSE, update.sigma=FALSE)
            setU(linit, init$u)
        }
        init <- linit
    } else if (is.function(init)) {
        init <- do.call(init,list(formula=formula, data=data, REML=TRUE, ...))
    } else if (is(init, "merMod") || is(init, "rlmerMod")) {
        ## check whether formula and data match with
        ## the ones in the provided in init
        if (is.null(icall <- getCall(init)))
            stop("Object 'init' should contain a 'call' component")
        if (!identical(as.character(lcall$formula), as.character(icall$formula)) |
            !identical(lcall$data, icall$data))
            warning("Arguments 'data' and 'formula' do not match with 'init': ",
                    "using model specification from 'init'")
    } else stop("Unsuitable init object, aborting.",
                "Expecting no, list (see ?rlmer), rlmerMod or merMod object")
    lobj <- as(init, "rlmerMod")
    lobj@call <- lcall

    ## give a warning if weights or offset are used
    if (any(lobj@resp$weights != 1))
        stop("Argument weights is unsave to use at the moment.")
    if (any(lobj@resp$offset != 0))
        warning("Argument offset is untested.")

    ## set arguments only relevant to rlmerMod
    ## convert rho argument to list if necessary
    convRho <- function(l)
        if (!is.list(l)) rep.int(list(l), length(lobj@dim)) else l
    lobj@rho.b <- rho.b <- convRho(rho.b)
    if (missing("rho.sigma.b")) {
        rho.sigma.b <- list()
        ## set default wExp.b
        wExp.b <- ifelse(lobj@dim == 1, 2, 1) ## if s==1 then 2 else 1
        ## check length of c.sigma.b
        for (bt in seq_along(lobj@blocks)) {
            ## set higher tuning constants by default
            
            rho.sigma.b[[bt]] <- switch(wExp.b[bt],
                                        rho.b[[bt]],
                                        psi2propII(rho.b[[bt]]),
                                        stop("only wExp = 1 and 2 are supported"))
        }
    }
    lobj@rho.sigma.b <- convRho(rho.sigma.b)
    if (!isTRUE(chk <- validObject(lobj))) stop(chk)
    lobj@rho.e <- rho.e
    if (missing("rho.sigma.e"))
        rho.sigma.e <- psi2propII(rho.e) ## prop II is the default
    lobj@rho.sigma.e <- rho.sigma.e
    if (method == "DAStau" & any(sapply(lobj@idx, nrow) > 2)) {
        warning("Method 'DAStau' does not support blocks of size larger than 2. ",
                "Falling back to method 'DASvar'.")
        method <- "DASvar"
    }
    lobj@method <- method
    if (substr(method, 1, 3) == "DAS") {
        lobj@pp <- as(lobj@pp, "rlmerPredD_DAS")
        lobj@pp$method <- method
    }
    lobj@pp$initRho(lobj)
    lobj@pp$initMatrices(lobj)
    lobj@pp$updateMatrices()


    if (!doFit) return(updateWeights(lobj))
    
    ## do not start with theta == 0
    if (any(theta(lobj)[lobj@lower == 0] == 0)) {
        if (verbose > 0)
            cat("Setting variance components from 0 to 1\n")
        theta0 <- theta(lobj)
        theta0[lobj@lower == 0 & theta0 == 0] <- 1
        setTheta(lobj, theta0, fit.effects = TRUE, update.sigma = method != "Opt")
    } else {
        ## set theta at least once
        setTheta(lobj, theta(lobj), fit.effects = FALSE)
    }

    if (verbose > 0) {
        cat("\nrlmer starting values:\n")
        cat("sigma, theta: ", lobj@pp$sigma, ", ", theta(lobj), "\n")
        cat("coef: ", lobj@pp$beta, "\n")
        if (verbose > 1)
            cat("b.s: ", b.s(lobj), "\n")
    }

    ## required for max.iter:
    r <- len(lobj, "theta")

    ## do fit: non diagonal case differently
    if (!isDiagonal(lobj@pp$U_b)) {
        if (method == "DASvar") {
            lobj <- rlmer.fit.DAS.nondiag(lobj, verbose, max.iter, rel.tol)
        } else if (method == "DAStau") {
            lobj <- rlmer.fit.DAS.nondiag(lobj, verbose, max.iter, rel.tol, method="DAStau")
        } else 
            stop("Non-diagonal case only supported by DAStau and DASvar")
    } else {
        lobj <- rlmer.fit.DAS(lobj, verbose, max.iter, rel.tol)
    }

    if (verbose > 0) {
        cat("sigma, theta: ", lobj@pp$sigma, ", ", theta(lobj), "\n")
        cat("coef: ", lobj@pp$beta, "\n")
        if (verbose > 1)
            cat("b.s: ", b.s(lobj), "\n")
    }

    return(updateWeights(lobj))
}

## DAS method
rlmer.fit.DAS.nondiag <- function(lobj, verbose, max.iter, rel.tol, method=lobj@method,
                                  checkFalseConvergence = TRUE) {
    if (!.isREML(lobj))
        stop("can only do REML when using averaged DAS-estimate for sigma")

    ## Prepare for DAStau
    if (method == "DAStau") {
        ## 4d int
        ## vectorize it!
        ghZ <- as.matrix(expand.grid(lobj@pp$ghz, lobj@pp$ghz, lobj@pp$ghz, lobj@pp$ghz))
        ghw <- apply(as.matrix(expand.grid(lobj@pp$ghw, lobj@pp$ghw, lobj@pp$ghw, lobj@pp$ghw)), 1, prod)
    }

    ## fit model using EM algorithm
    conv <- FALSE
    convBlks <- rep(FALSE, length(lobj@blocks))
    iter <- 0
    rel.tol <- sqrt(rel.tol)
    ## compute kappa
    kappas <- lobj@pp$kappa_b
    ## zero pattern for T matrix
    nzT <- crossprod(bdiag(lobj@blocks[lobj@ind])) == 0
    q <- lobj@pp$q
    ## false convergence indicator
    fc <- rep(FALSE, length(lobj@blocks))
    if (verbose > 0) {
        theta0 <- theta(lobj)
        if (verbose > 2) {
            coef0 <- lobj@pp$beta
            b.s0 <- b.s(lobj)
            sigma0 <- lobj@pp$sigma
        }
    }
    ## iterate
    while(!conv && (iter <- iter + 1) < max.iter) {
        if (verbose > 0) cat("---- Iteration", iter, " ----\n")
        thetatilde <- theta(lobj)
        sigma <- .sigma(lobj)

        ## get expected value of cov(\tvbs)
        q <- len(lobj, "b")
        T <- switch(method,
                    DASvar=lobj@pp$Tb(),
                    DAStau=calcTau.nondiag(lobj, ghZ, ghw, .S(lobj), kappas, max.iter,
                                           rel.tol = rel.tol, verbose = verbose),
                    stop("Non-diagonal case only implemented for DASvar"))
        ## compute robustness weights and add to t and bs
        T[nzT] <- 0
        ## symmetrize T to avoid non symmetric warning, then apply chol
        T <- symmpart(T)
        ## save to cache
        lobj@pp$setT(T)
        ## apply chol to non-zero part only
        idx <- !lobj@pp$zeroB
        ## stop if all are zero
        if (!any(idx)) break
        L <- t(chol(T[idx,idx]))
        T.bs <- numeric(q) ## set the others to zero
        T.bs[idx] <- forwardsolve(L, lobj@pp$b.s[idx])
        ## compute weights
        db <- .dk(lobj, sigma, FALSE, T.bs)[lobj@k]
        wbsEta <- wbsDelta <- numeric(0)
        for (type in seq_along(lobj@blocks)) {
            s <- lobj@dim[type]
            lidx <- as.vector(lobj@idx[[type]])
            if (s > 1) {
                ## for eta, we would actually need a little smaller
                ## tuning constants than for delta to get the same efficiency
                wbsEta <- c(wbsEta, lobj@rho.sigma.b[[type]]@wgt(db[lidx]))
                wbsDelta <- c(wbsDelta, (lobj@rho.sigma.b[[type]]@psi(db[lidx]) -
                                         lobj@rho.sigma.b[[type]]@psi(db[lidx] - s*kappas[type]))/s)
            } else {
                lw <- lobj@rho.sigma.b[[type]]@wgt(db[lidx])
                wbsEta <- c(wbsEta, lw)
                wbsDelta <- c(wbsDelta, lw*kappas[type]) ## adding kappa to wbsDelta in 1d case
            }
        }
        WbDelta <- Diagonal(x=wbsDelta)
        T <- WbDelta %*% T
        bs <- sqrt(wbsEta) * lobj@pp$b.s

        ## cycle block types
        for(type in seq_along(lobj@blocks)) {
            if (convBlks[type]) next
            bidx <- lobj@idx[[type]]
            if (verbose > 5) {
                cat("Tau for blocktype ", type, ":", as.vector(T[bidx[,1],bidx[,1]]), "\n")
            }
            ## catch dropped vc
            if (all(abs(bs[bidx]) < 1e-7)) {
                if (verbose > 1)
                    cat("Block", type, "dropped (all = 0), stopping iterations.\n")
                 Ubtilde <- lobj@blocks[[type]]
                pat <- Ubtilde != 0
                Lind <- Ubtilde[pat]
                thetatilde[Lind] <- 0
                convBlks[type] <- TRUE
                next
            }
            s <- nrow(bidx)
            K <- ncol(bidx)
            ## right hand side
            ## sum over blocks
            rhs <- matrix(0, s, s)
            for (k in 1:K) {
                ## add weights
                rhs <- rhs + T[bidx[,k],bidx[,k]]
            }
            rhs <- rhs / K
            ## left hand side
            lbs <- matrix(bs[bidx], ncol(bidx), nrow(bidx), byrow=TRUE)
            ## add left hand side
            lhs <- crossprod(lbs / sigma) / K
            if (verbose > 2) {
                cat("LHS:", as.vector(lhs), "\n")
                cat("RHS:", as.vector(rhs), "\n")
                cat("sum(abs(LHS - RHS)):", sum(abs(lhs - rhs)), "\n")
            }
            ## if (isTRUE(all.equal(rhs, lhs, check.attributes=FALSE, tolerance = rel.tol))) {
            diff <- abs(rhs - lhs)
            if (all(diff < rel.tol * max(diff, rel.tol))) {
                if (verbose > 1)
                    cat("Estimating equations satisfied for block", type,
                        ", stopping iterations.\n")
                convBlks[type] <- TRUE
                next
            }
            deltaT <- as(backsolve(lchol(rhs), lchol(lhs)), "sparseMatrix")
            if (verbose > 1) cat("deltaT:", deltaT@x, "\n")
            ## get old parameter estimates for this block
            Ubtilde <- lobj@blocks[[type]]
            pat <- Ubtilde != 0
            Lind <- Ubtilde[pat]
            diagLind <- diag(Ubtilde)
            Ubtilde[pat] <- thetatilde[Lind]
            ## update Ubtilde by deltaT
            thetatilde[Lind] <- tcrossprod(Ubtilde, deltaT)[pat]
            ## FIXME: check boundary conditions?
            ## check if varcomp is dropped
            if (all(thetatilde[diagLind] < 1e-7)) {
                thetatilde[Lind] <- 0
                convBlks[type] <- TRUE
                next
            }
            ## check if this block is converged
            diff <- abs(thetatilde[Lind] - theta(lobj)[Lind])
            if (verbose > 3)
                cat("criterion:", sum(diff), ">=",
                    rel.tol * max(diff, rel.tol), ":",
                    sum(diff) < rel.tol * max(diff, rel.tol), "\n")
            if (sum(diff) < rel.tol * max(diff, rel.tol)) {
                convBlks[type] <- TRUE
                ## check if estimating equations are satisfied
                if (checkFalseConvergence) {
                    if (verbose > 3)
                        cat("checking estimating equations:", sum(abs(lhs - rhs)),
                            ">", sqrt(rel.tol), ":", sum(abs(lhs - rhs)) > sqrt(rel.tol), "\n")
                    if (sum(abs(lhs - rhs)) > sqrt(rel.tol))
                        fc[type] <- TRUE
                }
                next
            }
        }
        ## set theta
        setTheta(lobj, thetatilde, fit.effects = TRUE,
                 update.sigma = FALSE)
        ## update sigma without refitting effects
        updateSigma(lobj, fit.effects = FALSE)
        if (verbose > 0) {
            cat("delta theta:", format(theta0 - thetatilde, nsmall=20, scientific=FALSE),
                "\n")
            theta0 <- thetatilde
            if (verbose > 1) {
                cat(sprintf("delta coef:  %.12f\n", sum(abs(coef0 - lobj@pp$beta))))
                cat(sprintf("delta u:     %.12f\n", sum(abs(b.s0 - b.s(lobj)))))
                cat(sprintf("delta sigma: %.12f\n", abs(sigma0 - lobj@pp$sigma)))
                coef0 <- lobj@pp$beta
                b.s0 <- b.s(lobj)
                sigma0 <- lobj@pp$sigma
                if (verbose > 2) {
                    cat("theta:", format(thetatilde, nsmall=20, scientific=FALSE), "\n")
                    cat("coef:   ", lobj@pp$beta,"\n")
                    cat("b.s:    ", b.s(lobj), "\n")
                    cat("sigmae: ", lobj@pp$sigma, "\n")
                }
            }
        }
        if (all(convBlks)) conv <- TRUE
    }

    optinfo <- list(optimizer = "rlmer.fit.DAS.nondiag",
                    conv = list(opt = 0),
                    feval = iter,
                    warnings = list(),
                    val = diff)
    
    if (iter == max.iter) {
        warning(wt <- "iterations did not converge, returning unconverged estimate.")
        optinfo$warnings <- list(wt)
        optinfo$conv$opt <- 1
    }
    if (any(fc)) {
        warning(wt <- "algorithm converged, but estimating equations are not satisfied.")
        optinfo$warnings <- c(optinfo$warnings, list(wt))
        optinfo$conv$opt <- 2
    }

    lobj@optinfo <- optinfo

    lobj
}

## DAS method
rlmer.fit.DAS <- function(lobj, verbose, max.iter, rel.tol) {
    if (!.isREML(lobj))
        stop("can only do REML when using averaged DAS-estimate for sigma")

    lupdateTheta <- switch(lobj@method,
                           DASvar=,
                           DAStau=updateThetaTau,
                           stop("method not supported by rlmer.fit.DAS:", lobj@method))

    ## fit
    converged <- FALSE
    theta0 <- theta(lobj)
    if (verbose > 1) {
        coef0 <- lobj@pp$beta
        b.s0 <- b.s(lobj)
        sigma0 <- lobj@pp$sigma
    }
    iter <- 0
    while (!converged && iter < max.iter) {
        iter <- iter + 1
        if (verbose > 0) cat("Iteration", iter, "\n")
        ## fit theta
        lupdateTheta(lobj, max.iter, rel.tol/10, verbose)
        theta1 <- theta(lobj)

        if (verbose > 0) {
            cat(sprintf("delta theta: %.12f\n", sum(abs(theta0 - theta1))))
            if (verbose > 1) {
                cat(sprintf("delta coef:  %.12f\n", sum(abs(coef0 - lobj@pp$beta))))
                cat(sprintf("delta u:     %.12f\n", sum(abs(b.s0 - b.s(lobj)))))
                cat(sprintf("delta sigma: %.12f\n", abs(sigma0 - lobj@pp$sigma)))
                coef0 <- lobj@pp$beta
                b.s0 <- b.s(lobj)
                sigma0 <- lobj@pp$sigma
                if (verbose > 2) {
                    cat("theta:  ", theta(lobj),"\n")
                    cat("coef:   ", lobj@pp$beta,"\n")
                    cat("b.s:    ", b.s(lobj), "\n")
                    cat("sigmae: ", lobj@pp$sigma, "\n")
                }
            }
        }

        ## all zero or change smaller than relative tolerance
        ## all zero: we can't get out of this anyway, so we have to stop.
        converged <- all(theta1 == 0) || sum(abs(theta0 - theta1)) <
            200*rel.tol*sum(abs(theta0))
        if (verbose > 1)
            cat(sprintf("Criterion: %.12f, %.12f", sum(abs(theta0 - theta1)),
                        sum(abs(theta0 - theta1)) / rel.tol / sum(abs(theta0)) / 200), "\n")
        theta0 <- theta1
    }

    optinfo <- list(optimizer = "rlmer.fit.DAS",
                    conv = list(opt = 0),
                    feval = iter,
                    warnings = list(),
                    val = diff)
    
    if (iter == max.iter) {
        warning(wt <- "iterations did not converge, returning unconverged estimate.")
        optinfo$warnings <- list(wt)
        optinfo$conv$opt <- 1
    }

    lobj@optinfo <- optinfo
    
    lobj
}
