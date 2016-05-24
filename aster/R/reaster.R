
reaster <- function(fixed, random, pred, fam, varvar, idvar, root,
    famlist = fam.default(), origin, data, effects, sigma, response)
    UseMethod("reaster")

reaster.default <- function(fixed, random, pred, fam, varvar, idvar, root,
    famlist = fam.default(), origin, data, effects, sigma, response)
{
    stopifnot(is.matrix(fixed))
    stopifnot(is.numeric(fixed))
    stopifnot(is.finite(fixed))

    if (! is.list(random))
        random <- list(random)
    for (irandom in seq(along = random)) {
        r <- random[[irandom]]
        if (! is.matrix(r))
            stop("random must be matrix or list of matrices")
        if (! is.numeric(r))
            stop("random must be numeric matrix or list of such")
        if (! all(is.finite(r)))
            stop("random must be all finite numeric matrix or list of such")
        if (nrow(r) != nrow(fixed))
            stop("fixed and random effect model matrices with different row dimensions")
    }

    stopifnot(is.vector(response))
    stopifnot(is.numeric(response))
    stopifnot(is.finite(response))
    if(length(response) != nrow(fixed))
        stop("length of response vector != row dimension of model matrices")

    stopifnot(is.vector(idvar))
    stopifnot(is.numeric(idvar))
    nind <- length(unique(idvar))

    stopifnot(is.factor(varvar))
    stopifnot(is.character(levels(varvar)))
    nnode <- nlevels(varvar)

    if (nind * nnode != length(response))
        stop("idvar or varvar wrong number of levels")

    varvarmat <- matrix(as.character(varvar), nind, nnode)
    idvarmat <- matrix(as.vector(idvar), nind, nnode)
    foo <- apply(varvarmat, 2, function(x) length(unique(x)))
    bar <- apply(idvarmat, 1, function(x) length(unique(x)))
    if (! (all(foo == 1) & all(bar == 1)))
        stop("idvar or varvar wrong structure")

    stopifnot(is.vector(pred))
    stopifnot(is.numeric(pred))
    stopifnot(pred %in% seq(0, nlevels(varvar)))

    stopifnot(is.vector(fam))
    stopifnot(is.numeric(fam))
    stopifnot(fam %in% seq(1, length(famlist)))
    famlist.blurfle <- sapply(famlist, as.character)
    famlist.problem <- grepl("negative.binomial", famlist.blurfle)
    if (any(famlist.problem[fam]))
        warning("negative binomial incompatible with random effects,",
            " see help(reaster)")

    stopifnot(is.vector(root))
    stopifnot(length(root) == length(response))
    stopifnot(is.numeric(root))
    stopifnot(is.finite(root))
    stopifnot(root > 0)

    if (! missing(origin)) {
        stopifnot(is.vector(origin))
        stopifnot(length(origin) == length(response))
        stopifnot(is.numeric(origin))
        stopifnot(is.finite(origin))
    }

    if (! missing(sigma)) {
        stopifnot(is.vector(sigma))
        stopifnot(is.numeric(sigma))
        stopifnot(is.finite(sigma))
        stopifnot(length(sigma) == length(as.list(random)))
    }

    ##### fit fixed effect aster model (even if effects supplied)
    ncoef <- ncol(fixed)
    modmat <- array(fixed, dim = c(nind, nnode, ncoef))
    dimnames(modmat) <- list(NULL, NULL, colnames(fixed))
    y <- matrix(response, nind, nnode)
    asterargs <- list(x = y, root = matrix(root, nind, nnode),
        pred = pred, fam = fam, modmat = modmat, famlist = famlist)
    if (! missing(origin))
        asterargs$origin <- matrix(origin, nind, nnode)
    aout <- try(do.call(aster, asterargs))
    if (inherits(aout, "try-error") || (! aout$converged)) {
        asterargs$method <- "nlm"
        aout <- try(do.call(aster, asterargs))
        if (inherits(aout, "try-error") || (! aout$converged))
            stop("cannot get convergence in fixed effects model")
    }

    ##### figure out what columns to drop from fixed effects matrix
    ikeep <- match(names(aout$coefficients), colnames(fixed))
    dropped <- colnames(fixed)[- ikeep]
    nfix.orig <- ncol(fixed)
    fixed <- fixed[ , ikeep]
    nfix <- ncol(fixed)
    idx.fixed <- seq(1, nfix)

    nrand <- sapply(random, ncol)
    blurfle <- names(random)
    if (is.null(blurfle)) {
        blurfle <- paste("R", seq(along = random), sep = "")
    } else {
        if (length(unique(blurfle)) != length(blurfle))
            stop("duplicate names for argument random")
    }
    splitter <- factor(rep(blurfle, times = nrand), levels = blurfle)
    idx.random <- nfix + seq(1, sum(nrand))

    if (! missing(effects)) {
        stopifnot(is.vector(effects))
        stopifnot(is.numeric(effects))
        stopifnot(is.finite(effects))
        if (length(dropped) != 0)
            if (length(effects) == nfix.orig + sum(nrand)) {
                warning("shortening effects to agree with columns dropped from fixed")
                effects.fixed <- effects[1:nfix.orig]
                effects.random <- effects[- (1:nfix.orig)]
                effects.fixed <- effects.fixed[ikeep]
                effects <- c(effects.fixed, effects.random)
            }
        if (length(effects) != ncol(fixed) + sum(sapply(random, ncol)))
            stop("length(effects) not sum of column dimensions of fixed and random effects model matrices")
    }

    if (missing(sigma) | missing(effects)) {
        sigma.start <- rep(1, length(nrand))
        eff.start <- c(aout$coefficients, rep(0, sum(nrand)))
        tout <- try(trust(objfun = penmlogl, eff.start,
            rinit = 1, rmax = 10,
            sigma = sigma.start, fixed = fixed, random = random,
            obj = aout, iterlim = 1000), silent = TRUE)
        if (inherits(tout, "try-error") || (! tout$converged))
            stop("cannot get convergence in penalized likelihood step")
        ### FIX ME ###
        ### need to try another optimization method when trust fails

        eff <- tout$argument * tout$scale
        eff.split <- split(eff[idx.random], splitter)
        sigma.start <- sapply(eff.split, function(x) sqrt(mean(x^2)))
    } else {
        sigma.start <- sigma
        eff.start <- effects
    }

    cache <- new.env(parent = emptyenv())
    oout <- suppressWarnings(try(optim(sigma.start, pickle, parm = eff.start,
        fixed = fixed, random = random, obj = aout, y = y, cache = cache),
        silent = TRUE))
    # this used to be
    #
    #     if (inherits(oout, "try-error") || oout$convergence != 0)
    #         stop("step 1 part 1 (optim Nelder-Mead with pickle) failed")
    #
    # but Marcus Warwell had an example where this optim with method
    # Nelder-Mead failed with out$convergence = 10 (indicates degeneracy
    # of the Nelder-Mead simplex, so says the documentation), and it is
    # very unclear what we are supposed to do with that.  Anyhoo, we
    # are only using this result as a starting point for more polishing,
    # so I hope this does not matter.

    if (inherits(oout, "try-error"))
        stop("step 1 part 1 (optim Nelder-Mead with pickle) failed")

    sigma.mle <- oout$par
    trout <- trust(objfun = penmlogl, cache$parm, rinit = 1, rmax = 10,
        sigma = sigma.mle, fixed = fixed, random = random, obj = aout, y = y)
    if (! trout$converged)
        stop("step 1 part 2 (trust with penmlogl) failed")
    parm.mle <- trout$argument
    zwz.mle <- makezwz(sigma.mle, parm = parm.mle,
        fixed = fixed, random = random, obj = aout, y = y)
    save.iter <- NULL
    value.mle <- NaN

    # now iterate pickle3 and trust
    foo <- function(alphaceesigma)
        pickle3(alphaceesigma, fixed = fixed, random = random,
            obj = aout, y = y, zwz = zwz.mle, deriv = 2)
    repeat {
        tout <- trust(foo, c(parm.mle, sigma.mle), rinit = 1, rmax = 10,
            iterlim = 1000)
        if (! tout$converged)
            stop("step 2 (trust with pickle3) failed")
        save.iter <- c(save.iter, tout$iterations)
        sigma.old <- as.vector(sigma.mle)
        sigma.mle <- tout$argument[nfix + sum(nrand) + seq(along = nrand)]
        parm.mle <- tout$argument[seq(1, nfix + sum(nrand))]
        zwz.mle <- makezwz(sigma.mle, parm.mle,
            fixed = fixed, random = random, obj = aout, y = y)
        value.mle <- tout$value
        if (isTRUE(all.equal(sigma.old, as.vector(sigma.mle)))) break
    }

    alpha.mle <- parm.mle[1:nfix]
    c.mle <- parm.mle[- (1:nfix)]

    # fix up negative sigma

    idx <- rep(seq(along = sigma.mle), times = nrand)
    for (k in seq(along = sigma.mle))
        if (sigma.mle[k] < 0) {
            sigma.mle[k] <- (- sigma.mle[k])
            eek <- (idx == k)
            c.mle[eek] <- (- c.mle[eek])
    }

    a.mle <- rep(sigma.mle, times = nrand)
    b.mle <- a.mle * c.mle

    if (! is.null(colnames(fixed)))
        names(alpha.mle) <- colnames(fixed)
    if (! is.null(colnames(random)))
        names(sigma.mle) <- names(random)
    if (all(! is.null(sapply(random, colnames)))) {
        names(b.mle) <- Reduce(c, Map(colnames, random))
        names(c.mle) <- names(b.mle)
    }

    if (missing(origin)) origin <- NULL
    result <- list(obj = aout, fixed = fixed, random = random,
        dropped = dropped, sigma = sigma.mle, nu = sigma.mle^2,
        c = c.mle, b = b.mle, alpha = alpha.mle, zwz = zwz.mle,
        response = response, origin = origin,
        iterations = save.iter, counts = oout$counts,
        deviance = 2.0 * value.mle)
    class(result) <- c("reaster", "asterOrReaster")
    return(result)
}

reaster.formula <- function(fixed, random, pred, fam, varvar, idvar, root,
    famlist = fam.default(), origin, data, effects, sigma, response)
{
    stopifnot(inherits(fixed, "formula"))
    if (! is.list(random))
        random <- list(random)
    if (!  all(sapply(random, function(x) inherits(x, "formula"))))
        stop("random must be formula or list of formulas")
    stopifnot(missing(data) || is.data.frame(data))

    save.fixed <- fixed
    save.random <- random

    oldopt <- options(na.action = na.fail)
    on.exit(options(oldopt))

    ##### stuff copied from glm.R and not understood #####
    ##### see also http://developer.r-project.org/model-fitting-functions.txt

    call <- match.call()
    if(missing(data))
        data <- environment(fixed)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("fixed", "data", "varvar", "idvar", "root"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fred <- as.list(mf)
    names(fred)[2] <- "formula"
    mf <- as.call(fred)
    mf <- eval.parent(mf)

    mt <- attr(mf, "terms")
    if (missing(response))
        response <- model.response(mf, "numeric")
    if (is.empty.model(mt)) {
        stop("empty model")
    } else {
        fixed <- model.matrix(mt, mf)
    }
    rownames(fixed) <- NULL

    varvar <- mf[["(varvar)"]]
    idvar <- mf[["(idvar)"]]
    root <- mf[["(root)"]]

    ##### end of stuff copied from glm.R and not understood #####

    for (irandom in seq(along = random)) {
        f <- random[[irandom]]
        foo <- lm(f, data = data, method = "model.frame")
        bar <- model.matrix(f, data = foo)
        rownames(bar) <- NULL
        random[[irandom]] <- bar
    }

    result <- NextMethod("reaster", response = response)

    result$formula <- list(fixed = save.fixed, random = save.random)
    result$call <- call
    class(result) <- c("reaster.formula", "reaster", "asterOrReaster")
    return(result)
}

summary.reaster <- function(object, standard.deviation = TRUE, ...)
{
    stopifnot(inherits(object, "reaster"))
    stopifnot(inherits(object$obj, "aster"))
    stopifnot(is.logical(standard.deviation))
    stopifnot(length(standard.deviation) == 1)

    alpha <- object$alpha
    sigma <- object$sigma
    nu <- object$nu
    cee <- object$c
    bee <- object$b
    zwz <- object$zwz
    fixed <- object$fixed
    random <- object$random
    if (is.matrix(random))
        random <- list(random)
    nfix <- ncol(fixed)
    nrand <- sapply(random, ncol)

    obj <- object$obj
    y <- matrix(object$response, nrow = nrow(obj$x), ncol = ncol(obj$x))
    origin <- object$origin

    if (is.null(origin)) {
        iz <- is.zero(c(alpha, bee, nu), fixed, random, obj, y, zwz = zwz)
    } else {
        iz <- is.zero(c(alpha, bee, nu), fixed, random, obj, y,
            origin = origin, zwz = zwz)
    }
    nu[iz] <- 0
    sigma[iz] <- 0
    izbee <- rep(iz, times = nrand)
    bee[izbee] <- 0

    has.se <- c(rep(TRUE, length(alpha)), ! iz)
    if (all(iz)) {
        se.alpha <- sqrt(diag(solve(obj$fisher)))
        se.bee <- rep(NA_real_, length(bee))
        se.sigma <- rep(NA_real_, length(sigma))
        se.nu <- rep(NA_real_, length(nu))
        # want to return "subfish" created in other part
        # not defining subfish here is bug that causes crash
        # when this case (all variance components zero) occurs
        subfish <- obj$fisher
    } else {
        subrandom <- random[! iz]
        subnu <- nu[! iz]
        subbee <- bee[! izbee]
        subzwz <- zwz[! izbee, , drop = FALSE]
        subzwz <- subzwz[ , ! izbee, drop = FALSE]
        if (is.null(origin)) {
            qout <- quickle(c(alpha, subnu), subbee, fixed, subrandom,
                obj, y, zwz = subzwz, deriv = 2)
        } else {
            qout <- quickle(c(alpha, subnu), subbee, fixed, subrandom,
                obj, y, origin = origin, zwz = subzwz, deriv = 2)
        }
        subfish <- qout$hessian
        eout <- eigen(subfish, symmetric = TRUE)
        goodfish <- min(eout$values) > 0
        if (goodfish) {
            # was buggy.  Have to use eigen to invert if use eigen to
            # check if invertible !!!!!!!!
            subfish.inv <- eout$vectors %*% diag(1 / eout$values) %*%
                t(eout$vectors)
            se.subparm <- sqrt(diag(subfish.inv))
        } else {
            warning(paste("estimated Fisher information matrix not positive",
               "definite, making all standard errors infinite"))
            se.subparm <- rep(Inf, nrow(subfish))
        }
        se.parm <- rep(NA_real_, length(has.se))
        se.parm[has.se] <- se.subparm
        se.alpha <- se.parm[seq(along = alpha)]
        se.nu <- se.parm[- seq(along = alpha)]
        se.sigma <- se.nu / (2 * sigma)
    }

    foo <- alpha
    foo <- cbind(foo, se.alpha)
    foo <- cbind(foo, foo[ , 1] / foo[ , 2])
    foo <- cbind(foo, 2 * pnorm(- abs(foo[ , 3])))
    rownames(foo) <- colnames(object$fixed)
    colnames(foo) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

    bar <- sigma
    bar <- cbind(bar, se.sigma)
    bar <- cbind(bar, bar[ , 1] / bar[ , 2])
    bar <- cbind(bar, pnorm(- abs(bar[ , 3])))
    blurfle <- names(object$random)
    if (is.null(blurfle))
        blurfle <- paste("sigma", seq(along = object$random), sep = "")
    rownames(bar) <- blurfle
    colnames(bar) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)/2")

    baz <- nu
    baz <- cbind(baz, se.nu)
    baz <- cbind(baz, baz[ , 1] / baz[ , 2])
    baz <- cbind(baz, pnorm(- abs(baz[ , 3])))
    blurfle <- names(object$random)
    if (is.null(blurfle))
        blurfle <- paste("nu", seq(along = object$random), sep = "")
    rownames(baz) <- blurfle
    colnames(baz) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)/2")

    return(structure(list(alpha = foo, sigma = bar, nu = baz,
        object = object, standard.deviation = standard.deviation,
        fisher = subfish), class = "summary.reaster"))
}

print.summary.reaster <-
    function (x, digits = max(3, getOption("digits") - 3),
        signif.stars = getOption("show.signif.stars"), ...)
{
    if (! is.null(x$object$call)) {
        cat("\nCall:\n")
        cat(paste(deparse(x$object$call), sep="\n", collapse="\n"),
            "\n\n", sep="")
    }

    cat("\nFixed Effects:\n")
    printCoefmat(x$alpha, digits = digits,
        signif.stars = signif.stars, na.print = "NA", ...)

    if (x$standard.deviation) {
        cat("\nSquare Roots of Variance Components (P-values are one-tailed):\n")
        printCoefmat(x$sigma, digits = digits,
            signif.stars = signif.stars, na.print = "NA", ...)
    } else {
        cat("\nVariance Components (P-values are one-tailed):\n")
        printCoefmat(x$nu, digits = digits,
            signif.stars = signif.stars, na.print = "NA", ...)
    }

    if (! is.null(x$bee)) {
        cat("\nRandom Effects:\n")
        printCoefmat(x$bee, digits = digits,
            signif.stars = signif.stars, na.print = "NA", ...)
    }

    return(invisible(x))
}

anova.reaster <- function(object, ...)
{
    dotargs <- list(...)
    if (length(dotargs) == 0)
        stop("need at least two objects of class \"reaster\"")
    allargs <- c(list(object), dotargs)
    if (! all(sapply(allargs, function(x) inherits(x, "reaster"))))
        stop("some arguments not of class \"reaster\"")
    return(anova.reasterlist(allargs))
}

anova.reasterlist <- function(object, ...)
{
    stopifnot(is.list(object))
    stopifnot(length(object) >= 2)
    if (! all(sapply(object, function(x) inherits(x, "reaster"))))
        stop("some components not of class \"reaster\"")

    nmodels <- length(object)

    # attempt to check that models are nested
    # all models have same origin
    ok <- TRUE
    for (i in 2:nmodels) {
        o1 <- object[[i - 1]]$origin
        o2 <- object[[i]]$origin
        ok <- ok & identical(o1, o2)
    }
    if (! ok) warning("not same origin for all models, probably not nested")

    # all fixed effect model matrices have column labels
    ok <- TRUE
    for (i in 1:nmodels) {
        mf <- object[[i]]$fixed
        ok <- ok & (! is.null(colnames(mf)))
    }
    if (! ok) warning("no colnames for some model matrix for fixed effects, cannot check nesting")

    # all fixed effect model matrices column labels are nested
    ok <- TRUE
    for (i in 2:nmodels) {
        mf1 <- object[[i - 1]]$fixed
        mf2 <- object[[i]]$fixed
        ok <- ok & all(colnames(mf1) %in% colnames(mf2))
    }
    if (! ok) warning("colnames for model matrices for fixed effects not nested, models probably not nested")

    # all random effect model matrices have column labels
    ok <- TRUE
    for (i in 1:nmodels) {
        mr <- object[[i]]$random
        cr <- lapply(mr, colnames)
        ok <- ok & (! any(is.null(cr)))
    }
    if (! ok) warning("no colnames for some model matrix for random effects, cannot check nesting")

    # all random effect model matrices column labels are nested
    ok <- TRUE
    for (i in 2:nmodels) {
        mr1 <- object[[i - 1]]$random
        mr2 <- object[[i]]$random
        cr1 <- lapply(mr1, colnames)
        cr2 <- lapply(mr2, colnames)
        for (j1 in seq(along = cr1)) {
            ok.too <- FALSE
            for (j2 in seq(along = cr2)) {
                ok.too <- ok.too | identical(cr1[[j1]], cr2[[j2]])
            }
            ok <- ok & ok.too
        }
    }
    if (! ok) warning("colnames for model matrices for random effects not nested, models probably not nested")

    resdf.fix <- sapply(object, function(x) nrow(x$fixed))
    resdf.rand <- sapply(object, function(x) length(x$random))
    if (any(resdf.rand > 1)) stop("don't yet know how to compare models differing by two or more variance components")
    resdf <- resdf.fix + resdf.rand
    resdev <- sapply(object, function(x) x$deviance)
    table <- data.frame(resdf, resdev, c(NA, diff(resdf.fix)),
        c(NA, diff(resdf.rand)), c(NA, diff(resdev)))
    fred <- function(x) {
        foo <- x$fixed
        if (is.null(foo)) return("(no formulas)")
        result <- as.character(deparse(foo))
        for (i in seq(along = x$random)) {
            foo <- x$random[[i]]
            if (is.null(foo)) return("(no formulas)")
            more.result <- as.character(deparse(foo))
            result <- paste(result, more.result, sep = ", ")
        }
        return(result)
    }
    variables <- sapply(object, fred)
    dimnames(table) <- list(1:nmodels, c("Model Df", "Model Dev", "Df Fix",
        "Df Rand", "Deviance"))
    title <- "Analysis of Deviance Table\n"
    topnote <- paste("Model ", format(1:nmodels),": ",
        variables, collapse = "\n", sep = "")
    mixers <- table[ , "Df Rand"] == 1
    pval <- double(nrow(table))
    pval[! mixers] <- pchisq(table[ , "Deviance"],
        table[ , "Df Fix"], lower.tail = FALSE)
    pval[mixers] <- (0.5) * pchisq(table[ , "Deviance"],
        table[ , "Df Fix"], lower.tail = FALSE) +
        pchisq(table[ , "Deviance"], table[ , "Df Fix"] + 1,
        lower.tail = FALSE)
    table <- cbind(table, "P-value" = pval)
    structure(table, heading = c(title, topnote),
              class = c("anova", "data.frame"))
}

