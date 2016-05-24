
aster <- function(x, ...)
    UseMethod("aster")

aster.default <- function(x, root, pred, fam, modmat, parm,
    type = c("unconditional", "conditional"), famlist = fam.default(),
    origin, origin.type = c("model.type", "unconditional", "conditional"),
    method = c("trust", "nlm", "CG", "L-BFGS-B"), fscale, maxiter = 1000,
    nowarn = TRUE, newton = TRUE, optout = FALSE, coef.names, ...)
{
    type <- match.arg(type)
    origin.type <- match.arg(origin.type)
    method <- match.arg(method)

    stopifnot(is.numeric(x))
    stopifnot(is.matrix(x))
    nind <- nrow(x)
    nnode <- ncol(x)
    stopifnot(is.numeric(root))
    stopifnot(is.matrix(root))
    stopifnot(identical(dim(x), dim(root)))
    stopifnot(length(pred) == ncol(x))
    stopifnot(all(pred == as.integer(pred)))
    stopifnot(all(pred < seq(along = pred)))
    stopifnot(length(fam) == ncol(x))
    stopifnot(all(fam == as.integer(fam)))
    stopifnot(all(is.element(fam, seq(along = famlist))))
    stopifnot(is.numeric(modmat))
    stopifnot(is.array(modmat))
    stopifnot(length(dim(modmat)) == 3)
    stopifnot(identical(dim(modmat)[1:2], dim(x)))

    if (missing(parm)) {
        parm <- rep(0, dim(modmat)[3])
    } else {
        parm <- as.double(parm)
        if (length(parm) != dim(modmat)[3])
            stop("parm wrong length, not dimension 3 of modmat")
    }

    setfam(famlist)

    if (missing(origin)) {
        origin <- .C("aster_default_origin",
            nind = as.integer(nind),
            nnode = as.integer(nnode),
            fam = as.integer(fam),
            theta = matrix(as.double(0), nind, nnode),
            PACKAGE = "aster")$theta
        if (type == "unconditional")
            origin <- .C("aster_theta2phi",
                nind = as.integer(nind),
                nnode = as.integer(nnode),
                pred = as.integer(pred),
                fam = as.integer(fam),
                theta = origin,
                phi = matrix(as.double(0), nind, nnode),
                PACKAGE = "aster")$phi
    } else {
        stopifnot(is.numeric(origin))
        storage.mode(origin) <- "double"
        stopifnot(length(origin) == nind * nnode)
        stopifnot(all(is.finite(origin)))
        if (is.matrix(origin)) {
            stopifnot(nrow(origin) == nind)
            stopifnot(ncol(origin) == nnode)
        } else {
            origin <- matrix(origin, nind, nnode)
        }
        if (origin.type == "model.type")
            origin.type <- type
        if (type == "unconditional" && origin.type == "conditional")
            origin <- .C("aster_theta2phi",
                nind = as.integer(nind),
                nnode = as.integer(nnode),
                pred = as.integer(pred),
                fam = as.integer(fam),
                theta = origin,
                phi = matrix(as.double(0), nind, nnode),
                PACKAGE = "aster")$phi
        if (type == "conditional" && origin.type == "unconditional")
            origin <- .C("aster_phi2theta",
                nind = as.integer(nind),
                nnode = as.integer(nnode),
                pred = as.integer(pred),
                fam = as.integer(fam),
                phi = origin,
                theta = matrix(as.double(0), nind, nnode),
                PACKAGE = "aster")$theta
    }

    ##### try starting parm and origin
    mout <- mloglhelper(parm, pred, fam, x, root, modmat, origin,
                deriv = 0, type = type)
    if (! is.finite(mout$value))
        stop("initial \"origin\", \"modmat\", and \"parm\" invalid")

    mtry <- matrix(as.numeric(modmat), nrow = nind * nnode)
    qtry <- qr(mtry)
    if (qtry$rank < ncol(mtry)) {
        inies <- qtry$pivot[seq(1, qtry$rank)]
        cnames <- dimnames(modmat)[[3]]
        outies <- cnames[- inies]
        modmat <- modmat[ , , inies, drop = FALSE]
        parm <- parm[inies]
    } else {
        outies <- character(0)
    }

    ncoef <- dim(modmat)[3]
    if (missing(fscale))
        fscale <- nind
    if (missing(coef.names))
        coef.names <- dimnames(modmat)[[3]]

    if (method == "trust") {
        objfun <- function(beta) {
            mloglhelper(beta, pred, fam, x, root, modmat, origin,
                deriv = 2, type = type)
        }
        otherargs <- list(...)
        rinit <- otherargs$rinit
        rmax <- otherargs$rmax
        if (is.null(rinit)) rinit <- 1
        if (is.null(rmax)) rmax <- 10
        if (nowarn) {
            suppressWarnings(oout <- trust(objfun, parm, rinit, rmax,
                iterlim = maxiter, ...))
        } else {
            oout <- trust(objfun, parm, rinit, rmax,
                iterlim = maxiter, ...)
        }
        aout <- list(coefficients = oout$argument,
            iter = oout$iterations, converged = oout$converged)
    }
    if (method == "nlm") {
        objfun <- function(beta) {
            out <- mloglhelper(beta, pred, fam, x, root, modmat, origin,
                deriv = 1, type = type)
            result <- out$value
            attr(result, "gradient") <- out$gradient
            return(result)
        }
        if (nowarn) {
            suppressWarnings(oout <- nlm(objfun, parm, fscale = fscale,
                iterlim = maxiter, ...))
        } else {
            oout <- nlm(objfun, parm, fscale = fscale,
                iterlim = maxiter, ...)
        }
        aout <- list(coefficients = oout$estimate,
            iter = oout$iterations, code = oout$code,
            converged = oout$code <= 2)
    }
    if (method == "CG" || method == "L-BFGS-B") {
        objfun <- function(beta)
            mloglhelper(beta, pred, fam, x, root, modmat, origin,
                deriv = 0, type = type)$value
        grdfun <- function(beta)
            mloglhelper(beta, pred, fam, x, root, modmat, origin,
                deriv = 1, type = type)$gradient
        have.control <- is.element("control", names(list(...)))
        if (have.control) {
            stopifnot(is.list(control))
            control$fnscale <- fscale
            control$maxit <- maxiter
        } else {
            control <- list(fnscale = fscale, maxit = maxiter)
        }
        if (nowarn) {
            suppressWarnings(oout <- optim(parm, objfun, grdfun,
                method = method, control = control, ...))
        } else {
            oout <- optim(parm, objfun, grdfun,
                method = method, control = control, ...)
        }
        aout <- list(coefficients = oout$par,
            iter = oout$counts, code = oout$convergence,
            converged = oout$convergence == 0)
        if (method == "L-BFGS-B")
            aout$message <- oout$message
    }
    if (optout)
        aout$optout <- oout
    mout <- mloglhelper(aout$coefficients, pred, fam, x, root, modmat,
        origin, deriv = 2, type = type)
    if (newton && method != "trust") {
        qux <- qr(mout$hessian)
        if (qux$rank < dim(mout$hessian)[1]) {
            warning("rank deficient hessian, Newton step skipped")
        } else {
            beta <- aout$coefficients - solve(qux, mout$gradient)
            mout <- mloglhelper(beta, pred, fam, x, root, modmat, origin,
                deriv = 2, type = type)
            aout$coefficients <- beta
        }
    }
    aout$deviance <- 2 * mout$value
    aout$gradient <- mout$gradient
    aout$hessian <- mout$hessian
    aout$newton <- newton
    aout$rank <- qr(mout$hessian)$rank
    aout$x <- x
    aout$root <- root
    aout$pred <- pred
    aout$fam <- fam
    aout$modmat <- modmat
    aout$type <- type
    aout$famlist <- famlist
    names(aout$coefficients) <- coef.names
    if (type == "conditional") {
        fout <- .C("aster_fish_cond",
            nind = as.integer(nind),
            nnode = as.integer(nnode),
            ncoef = as.integer(ncoef),
            pred = as.integer(pred),
            fam = as.integer(fam),
            beta = as.double(aout$coefficients),
            root = as.double(root),
            x = as.double(x),
            modmat = as.double(modmat),
            fish = matrix(as.double(0), ncoef, ncoef),
            PACKAGE = "aster")
        aout[["fisher"]] <- fout$fish
    } else {
        aout[["fisher"]] <- mout$hessian
    }
    class(aout) <- c("aster", "asterOrReaster")
    if (! aout$converged)
        warning("Algorithm did not converge")
    if (length(outies) > 0)
        aout$dropped <- outies
    aout$origin <- origin
    clearfam()
    return(aout)
}

aster.formula <- function(formula, pred, fam, varvar, idvar, root,
    data, parm, type = c("unconditional", "conditional"),
    famlist = fam.default(),
    origin, origin.type = c("model.type", "unconditional", "conditional"),
    method = c("trust", "nlm", "CG", "L-BFGS-B"), fscale, maxiter = 1000,
    nowarn = TRUE, newton = TRUE, optout = FALSE, ...)
{
    type <- match.arg(type)
    origin.type <- match.arg(origin.type)
    method <- match.arg(method)

    oldopt <- options(na.action = na.fail)
    on.exit(options(oldopt))

    ##### stuff copied from glm.R and not understood #####
    ##### see also http://developer.r-project.org/model-fitting-functions.txt

    call <- match.call()
    if(missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "varvar", "idvar", "root"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)

    mt <- attr(mf, "terms")
    x <- model.response(mf, "numeric")
    if (is.empty.model(mt)) {
        stop("empty model")
    } else {
        modmat <- model.matrix(mt, mf)
    }

    varvar <- mf[["(varvar)"]]
    idvar <- mf[["(idvar)"]]
    root <- mf[["(root)"]]

    ##### end of stuff copied from glm.R and not understood #####

    nind <- length(unique(idvar))
    nnode <- length(unique(varvar))
    if (nind * nnode != length(varvar))
        stop("nrow(data) not nind * nnode")
    varvarmat <- matrix(as.vector(varvar), nind, nnode)
    idvarmat <- matrix(as.vector(idvar), nind, nnode)
    foo <- apply(varvarmat, 2, function(x) length(unique(x)))
    bar <- apply(idvarmat, 1, function(x) length(unique(x)))
    if (! (all(foo == 1) & all(bar == 1)))
        stop("data not nind by nnode matrix with rows individuals and columns variables")
    varlab <- varvarmat[1, ]
    idlab <- idvarmat[ , 1]
    if (all(idlab == seq(along = idlab)))
        idlab <- NULL

    if (! is.numeric(x))
        stop("response not numeric")
    if (length(x) != nind * nnode)
        stop("response not nind by nnode matrix with rows individuals and columns variables")
    x <- matrix(x, nind, nnode)
    dimnames(x) <- list(idlab, varlab)

    if (! is.numeric(root))
        stop("root not numeric")
    if (length(root) != nind * nnode)
        stop("root not nind by nnode matrix with rows individuals and columns variables")
    root <- matrix(root, nind, nnode)
    dimnames(root) <- list(idlab, varlab)

    if (! is.numeric(modmat))
        stop("model matrix not numeric")
    if (! is.matrix(modmat))
        stop("model matrix not matrix")
    if (nrow(modmat) != nind * nnode)
        stop("nrow of model matrix not nind * nnode")
    ncoef <- ncol(modmat)
    coeflab <- dimnames(modmat)[[2]]
    modmat <- array(as.numeric(modmat), c(nind, nnode, ncoef))
    dimnames(modmat) <- list(idlab, varlab, coeflab)

    if (missing(parm))
        parm <- rep(0, ncoef)
    if (missing(fscale))
        fscale <- nind

    out <- aster(x, root, pred, fam, modmat, parm, type, famlist, origin,
        origin.type, method, fscale, maxiter, nowarn, newton, optout, ...)

    class(out) <- c("aster.formula", "aster", "asterOrReaster")
    out$call <- call
    out$formula <- formula
    out$terms <- mt
    out$data <- data
    out$xlevels <- .getXlevels(mt, mf)
    return(out)
}

summary.aster <- function(object, info = c("expected", "observed"),
    info.tol = sqrt(.Machine$double.eps), show.graph = FALSE, ...)
{
    info <- match.arg(info)
    if (info == "expected")
        infomat <- object$fisher
    else
        infomat <- object$hessian

    if (! object$converged)
        stop("aster model fit not converged")

    fred <- eigen(infomat, symmetric = TRUE)
    sally <- fred$values < max(fred$values) * info.tol
    if (any(sally)) {
        cat("apparent null eigenvectors of information matrix\n")
        cat("directions of recession or constancy of log likelihood\n")
        print(zapsmall(fred$vectors[ , sally]))
        stop("cannot compute standard errors")
    }

    foo <- object$coefficients
    foo <- cbind(foo, sqrt(diag(solve(infomat))))
    foo <- cbind(foo, foo[ , 1] / foo[ , 2])
    foo <- cbind(foo, 2 * pnorm(- abs(foo[ , 3])))
    dimnames(foo) <- list(names(object$coefficients),
        c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    object$coefficients <- foo

    if (show.graph) {
        foo <- dimnames(object$x)[[2]]
        foo <- cbind(foo, c("root", foo)[object$pred + 1])
        bar <- object$famlist[object$fam]
        foo <- cbind(foo, sapply(bar, as.character))
        dimnames(foo) <- list(rep("", nrow(foo)),
            c("variable", "predecessor", "family"))
        object$graph <- foo
    }

    class(object) <- "summary.aster"
    return(object)
}

print.summary.aster <-
    function (x, digits = max(3, getOption("digits") - 3),
        signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")
    if (! is.null(x$graph)) {
        cat("\nGraphical Model:\n")
        print(x$graph, quote = FALSE)
        cat("\n")
    }
    printCoefmat(x$coefficients, digits = digits,
        signif.stars = signif.stars, na.print = "NA", ...)
    if (! is.null(x$dropped)) {
        cat("\nOriginal predictor variables dropped (aliased)\n")
        for (foo in x$dropped)
            cat("    ", foo, "\n")
    }
    return(invisible(x))
}

anova.aster <- function(object, ...)
{
    dotargs <- list(...)
    if (length(dotargs) == 0)
        stop("need at least two objects of class \"aster\"")
    if (! all(sapply(dotargs, function(x) inherits(x, "aster"))))
        stop("some arguments not of class \"aster\"")
    return(anova.asterlist(c(list(object), dotargs)))
}

anova.asterlist <- function(object, ...)
{
    stopifnot(is.list(object))
    stopifnot(length(object) >= 2)
    if (! all(sapply(object, function(x) inherits(x, "aster"))))
        stop("some components not of class \"aster\"")

    nmodels <- length(object)
    resdf  <- as.numeric(lapply(object, function(x) length(x$coefficients)))
    resdev <- as.numeric(lapply(object, function(x) x$deviance))
    table <- data.frame(resdf, resdev, c(NA, diff(resdf)),
        c(NA, - diff(resdev)))
    variables <- sapply(object, function(x) ifelse(is.null(x$formula),
        "(no formula)", as.character(deparse(x$formula))))
    dimnames(table) <- list(1:nmodels, c("Model Df", "Model Dev", "Df",
        "Deviance"))
    title <- "Analysis of Deviance Table\n"
    topnote <- paste("Model ", format(1:nmodels),": ",
        variables, collapse = "\n", sep = "")
    table <- cbind(table, "P(>|Chi|)" = pchisq(table[ , "Deviance"],
        table[ , "Df"], lower.tail = FALSE))
    structure(table, heading = c(title, topnote),
              class = c("anova", "data.frame"))
}

