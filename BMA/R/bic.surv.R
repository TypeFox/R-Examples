"bic.surv" <-
function (x, ...) 
UseMethod("bic.surv")
"bic.surv.data.frame" <-
function (x, surv.t, cens, strict = FALSE, OR = 20, maxCol = 30, 
    prior.param = c(rep(0.5, ncol(x))), OR.fix = 2, nbest = 150, 
    factor.type = TRUE, factor.prior.adjust = FALSE, call = NULL, ...) 
{
    leaps.bs <- function(info, coef, names.arg, nbest = nbest) {
        names.arg <- names.arg
        if (is.null(names.arg)) 
            names.arg <- c(as.character(1:9), LETTERS, letters)[1:ncol(info)]
        if (length(names.arg) < ncol(info)) 
            stop("Too few names")
        bIb <- coef %*% info %*% coef
        kx <- ncol(info)
        maxreg <- nbest * kx
        if (kx < 3) 
            stop("Too few independent variables")
        imeth <- 1
        df <- kx + 1
        Ib <- info %*% coef
        rr <- cbind(info, Ib)
        rr <- rbind(rr, c(Ib, bIb))
        it <- 0
        n.cols <- kx + 1
        nv <- kx + 1
        nf <- 0
        no <- 1e+05
        ib <- 1
        mb <- nbest
        nd <- n.cols
        nc <- 4 * n.cols
        rt <- matrix(rep(0, times = nd * nc), ncol = nc)
        rt[, 1:n.cols] <- rr
        iw <- c(1:(kx + 1), rep(0, times = 4 * nd))
        nw <- length(iw)
        rw <- rep(0, times = 2 * mb * kx + 7 * nd)
        nr <- length(rw)
        t1 <- 2
        s2 <- -1
        ne <- 0
        iv <- 0
        nret <- mb * kx
        Subss <- rep(0, times = nret)
        RSS <- Subss
        ans <- .Fortran("fwleaps", as.integer(nv), as.integer(it), 
            as.integer(kx), as.integer(nf), as.integer(no), as.integer(1), 
            as.double(2), as.integer(mb), as.double(rt), as.integer(nd), 
            as.integer(nc), as.integer(iw), as.integer(nw), as.double(rw), 
            as.integer(nr), as.double(t1), as.double(s2), as.integer(ne), 
            as.integer(iv), as.double(Subss), as.double(RSS), 
            as.integer(nret), PACKAGE = "BMA")
        regid <- ans[[21]]/2
        r2 <- ans[[20]]
        nreg <- sum(regid > 0)
        regid <- regid[1:nreg]
        r2 <- r2[1:nreg]
        which <- matrix(TRUE, nreg, kx)
        z <- regid
        which <- matrix(as.logical((rep.int(z, kx)%/%rep.int(2^((kx - 
            1):0), rep.int(length(z), kx)))%%2), byrow = FALSE, ncol = kx)
        size <- which %*% rep(1, kx)
        label <- character(nreg)
        sep <- if (all(nchar(names.arg) == 1)) 
            ""
        else ","
        for (i in 1:nreg) label[i] <- paste(names.arg[which[i, 
            ]], collapse = sep)
        ans <- list(r2 = r2, size = size, label = label, which = which)
        return(ans)
    }
    factor.names <- function(x) {
        out <- list()
        for (i in 1:ncol(x)) if (is.factor(x[, i])) 
            out[[i]] <- levels(x[, i])
        else out <- c(out, list(NULL))
        attributes(out)$names <- names(x)
        return(out)
    }
    dropcols <- function(x, surv.t, cens, maxCols = 30) {
        vnames <- attributes(x)$names
        nvar <- length(vnames)
        isfac <- rep(FALSE, times = nvar)
        for (i in 1:nvar) isfac[i] <- is.factor(x[, i])
        nlevels <- rep(NA, times = nvar)
        for (i in 1:nvar) if (isfac[i]) 
            nlevels[i] <- length(levels(x[, i]))
        any.dropped <- FALSE
        mm <- model.matrix(terms.formula(~., data = x), data = x)
        designx <- attributes(mm)$assign
        n.designx <- length(designx)
        designx.levels <- rep(1, times = n.designx)
        for (i in 2:n.designx) if (isfac[designx[i]]) 
            designx.levels[i] <- sum(designx[1:i] == designx[i]) + 
                1
        x.coxph <- data.frame(mm[, -1], surv.t = surv.t, cens = cens)
        cox.out <- coxph(Surv(surv.t, cens) ~ ., data = x.coxph, 
            method = "breslow", iter.max = 30)
        while (length(cox.out$coefficients) > maxCol) {
            any.dropped <- TRUE
            dropcox <- drop1(cox.out, test = "Chisq")
            dropped <- which.max(dropcox$"Pr(Chi)"[-1]) + 1
            x.coxph <- x.coxph[, -(dropped - 1)]
            designx.levels <- designx.levels[-dropped]
            designx <- designx[-dropped]
            cox.out <- coxph(Surv(surv.t, cens) ~ ., data = x.coxph, 
                method = "breslow", iter.max = 30)
        }
        remaining.vars <- unique(designx[-1])
        new.nvar <- length(remaining.vars)
        dropped.vars <- vnames[-remaining.vars]
        dropped.levels <- NULL
        ncol.cox <- ncol(x.coxph) - 2
        x.coxph <- x.coxph[-((ncol.cox + 1):(ncol.cox + 2))]
        xx <- data.frame(matrix(rep(NA, times = new.nvar * nrow(x.coxph)), 
            ncol = new.nvar))
        new.names = rep(NA, times = new.nvar)
        for (i in 1:new.nvar) {
            cvar <- remaining.vars[i]
            lvls <- designx.levels[cvar == designx]
            if (isfac[cvar]) {
                if (length(lvls) != length(levels(x[, cvar]))) {
                  newvar <- (as.matrix(x.coxph[, cvar == designx[-1]]) %*% 
                    cbind(lvls - 1)) + 1
                  xx[, i] <- factor(levels(x[, cvar])[newvar])
                  new.names[i] <- vnames[cvar]
                  removed.levels <- levels(x[, cvar])[-c(1, lvls)]
                  dropped.levels <- c(dropped.levels, paste(vnames[cvar], 
                    "_", removed.levels, sep = ""))
                }
                else {
                  xx[, i] <- factor(x[, cvar])
                  new.names[i] <- vnames[cvar]
                }
            }
            else {
                xx[, i] <- x[, cvar]
                new.names[i] <- vnames[cvar]
            }
        }
        dropped <- c(dropped.vars, dropped.levels)
        return(list(mm = xx, any.dropped = any.dropped, dropped = dropped, 
            var.names = new.names, remaining.vars = remaining.vars))
    }
    if (is.null(call)) 
        cl <- match.call()
    else cl <- call
    options(contrasts = c("contr.treatment", "contr.treatment"))
    x <- data.frame(x)
    prior.weight.denom <- 0.5^ncol(x)
    x.omit <- na.omit(x)
    used <- match(row.names(x), row.names(x.omit))
    omitted <- seq(nrow(x))[is.na(used)]
    if (length(omitted) > 0) {
        x <- x.omit
        surv.t <- surv.t[-omitted]
        cens <- cens[-omitted]
        warning(paste("There were ", length(omitted), "records deleted due to NA's"))
    }
    leaps.x <- x
    fac.levels <- rep(1, times = ncol(x))
    fn <- factor.names(x)
    factors <- !all(unlist(lapply(fn, is.null)))
    if (factors) {
        cdf <- cbind.data.frame(y = surv.t, x)
        mm <- model.matrix(formula(cdf), data = cdf)[, -1, drop = FALSE]
        mmm <- data.frame(matrix(mm, nrow = nrow(mm), byrow = FALSE))
        names(mmm) <- dimnames(mm)[[2]]
        output.names <- names(mmm)
        x.coxph <- data.frame(surv.t = surv.t, cens = cens, x)
        cox.out <- coxph(Surv(surv.t, cens) ~ ., data = x.coxph, 
            method = "breslow")
        cox.assign <- cox.out$assign
        for (tempi in length(cox.assign):1) {
            cox.assign[[tempi + 1]] <- cox.assign[[tempi]]
            names(cox.assign)[tempi + 1] <- names(cox.assign)[tempi]
        }
        cox.assign[[1]] <- numeric(0)
        names(cox.assign)[1] <- "(Intercept)"
        fac.levels <- unlist(lapply(cox.assign, length)[-1])
        if (factor.type) {
            for (i in 1:length(names(x))) {
                if (!is.null(fn[[i]])) {
                  nx <- names(x)[i]
                  coefs <- cox.out$coef[cox.out$assign[[i]]]
                  old.vals <- x[, i]
                  new.vals <- c(0, coefs)
                  new.vec <- as.vector(new.vals[match(old.vals, 
                    fn[[i]])])
                  leaps.x[, nx] <- new.vec
                }
            }
        }
        else {
            fac.levels <- rep(1, times = ncol(mmm))
            new.prior <- NULL
            for (i in 1:length(names(x))) {
                addprior <- prior.param[i]
                if (!is.null(fn[[i]])) {
                  k <- length(fn[[i]])
                  if (factor.prior.adjust) 
                    addprior <- rep(1 - (1 - prior.param[i])^(1/(k - 
                      1)), k - 1)
                  else addprior <- rep(prior.param[i], k - 1)
                }
                new.prior <- c(new.prior, addprior)
            }
            prior.param <- new.prior
            x <- leaps.x <- mmm
        }
    }


    xx <- data.frame()
    xx <- dropcols(leaps.x, surv.t, cens, maxCol)
    var.names <- xx$var.names
    remaining <- xx$remaining.vars
    leaps.x <- xx$mm
    reduced <- xx$any.dropped
    dropped <- NULL
    if (reduced) 
        dropped <- xx$dropped
    nvar <- length(x[1, ])
    x <- x[, remaining, drop = FALSE]
    x <- data.frame(x)
    fac.levels <- fac.levels[remaining]
    output.names <- list()
    for (i in 1:length(var.names)) {
        if (is.factor(x[, i])) 
            output.names[[i]] <- levels(x[, i])
        else output.names[[i]] <- NA
    }
    xnames <- names(x)
    names(leaps.x) <- var.names
    x.coxph <- data.frame(surv.t = surv.t, cens = cens, leaps.x)
    cox.out <- coxph(Surv(surv.t, cens) ~ ., data = x.coxph, 
        method = "breslow")
    x.coxph.fac <- data.frame(surv.t = surv.t, cens = cens, x)
    cox.assign <- coxph(Surv(surv.t, cens) ~ ., data = x.coxph.fac, 
        method = "breslow")$assign
    for (tempi in length(cox.assign):1) {
        cox.assign[[tempi + 1]] <- cox.assign[[tempi]]
        names(cox.assign)[tempi + 1] <- names(cox.assign)[tempi]
    }
    cox.assign[[1]] <- numeric(0)
    names(cox.assign)[1] <- "(Intercept)"
    n <- sum(cens)
    if (ncol(leaps.x) >= 3) {
        coef <- cox.out$coef
        info <- solve(cox.out$var)
        a <- leaps.bs(info, coef, names.arg = xnames, nbest = nbest)
        a$r2 <- c(0, a$r2)
        a$size <- c(0, a$size)
        a$label <- c("NULL", a$label)
        a$which <- rbind(rep(FALSE, ncol(x)), a$which)
        nmod <- length(a$size)
        prior.mat <- matrix(rep(prior.param, nmod), nmod, ncol(x), 
            byrow = TRUE)
        prior <- apply(a$which * prior.mat + (!a$which) * (1 - 
            prior.mat), 1, prod)
        bIb <- as.numeric(coef %*% info %*% coef)
        lrt <- bIb - (a$r2 * bIb)
        bic <- lrt + (a$size) * log(n) - 2 * log(prior)
        occam <- bic - min(bic) < 2 * OR.fix * log(OR)
        size <- a$size[occam]
        label <- a$label[occam]
        which <- a$which[occam, , drop = FALSE]
        bic <- bic[occam]
        prior <- prior[occam]
    }
    else {
        nmod <- switch(ncol(x), 2, 4)
        bic <- label <- rep(0, nmod)
        which <- matrix(c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE), nmod, nmod/2)
        size <- c(0, 1, 1, 2)[1:nmod]
        sep <- ","
        prior.mat <- matrix(rep(prior.param, nmod), nmod, ncol(x), 
            byrow = TRUE)
        prior <- apply(which * prior.mat + (!which) * (1 - prior.mat), 
            1, prod)
        for (k in 1:nmod) {
            if (k == 1) 
                label[k] <- "NULL"
            else label[k] <- paste(names(leaps.x)[which[k, ]], 
                collapse = sep)
        }
    }
    model.fits <- as.list(rep(0, length(label)))
    loglik <- rep(0, length(label))
    size <- rep(0, length(label))
    loglik.null <- cox.out$loglik[1]
    for (k in (1:length(label))) {
        if (sum(which[k, ]) != 0) {
            x.coxph <- data.frame(x[, which[k, ], drop = FALSE], 
                surv.t = surv.t, cens = cens)
            cox.out <- coxph(Surv(surv.t, cens) ~ ., data = x.coxph, 
                iter.max = 30, method = "breslow")
            loglik[k] <- cox.out$loglik[2]
            size[k] <- length(cox.out$coef)
            model.fits[[k]] <- matrix(rep(0, 2 * length(cox.out$coef)), 
                ncol = 2)
            model.fits[[k]][, 1] <- cox.out$coef
            model.fits[[k]][, 2] <- sqrt(diag(cox.out$var))
        }
        else {
            loglik[k] <- loglik.null
        }
    }
    bic <- size * log(n) - 2 * (loglik - loglik.null) - 2 * log(prior)
    occam <- bic - min(bic) < 2 * log(OR)
    size <- size[occam]
    label <- label[occam]
    which <- which[occam, , drop = FALSE]
    bic <- bic[occam]
    prior <- prior[occam]
    model.fits <- model.fits[occam]
    postprob <- (exp(-0.5 * (bic - min(bic))))/sum(exp(-0.5 * 
        (bic - min(bic))))
    order.bic <- order(bic, size, label)
    size <- size[order.bic]
    label <- label[order.bic]
    which <- which[order.bic, , drop = FALSE]
    bic <- bic[order.bic]
    prior <- prior[order.bic]
    postprob <- postprob[order.bic]
    model.fits <- model.fits[order.bic]
    nmod <- length(size)
    if (strict & (nmod != 1)) {
        occam <- rep(TRUE, nmod)
        for (k in (2:nmod)) {
            for (j in (1:(k - 1))) {
                which.diff <- which[k, ] - which[j, ]
                if (all(which.diff >= 0)) 
                  occam[k] <- FALSE
            }
        }
        size <- size[occam]
        label <- label[occam]
        which <- which[occam, , drop = FALSE]
        bic <- bic[occam]
        prior <- prior[occam]
        postprob <- postprob[occam]/(sum(postprob[occam]))
        model.fits <- model.fits[occam]
    }
    bic <- bic + 2 * log(prior)
    probne0 <- round(100 * t(which) %*% as.matrix(postprob), 
        1)
    nmod <- length(bic)
    nvar <- max(unlist(cox.assign))
    Ebi <- rep(0, nvar)
    SDbi <- rep(0, nvar)
    EbiMk <- matrix(rep(0, nmod * nvar), nrow = nmod)
    sebiMk <- matrix(rep(0, nmod * nvar), nrow = nmod)
    for (i in (1:ncol(x))) {
        whereisit <- cox.assign[[i + 1]]
        if (any(which[, i])) 
            for (k in (1:nmod)) if (which[k, i] == TRUE) {
                spot <- sum(which[k, (1:i)])
                posMk <- (c(0, cumsum(fac.levels[which[k, ]])) + 
                  1)[spot]
                posMk <- posMk:(posMk + fac.levels[i] - 1)
                EbiMk[k, whereisit] <- model.fits[[k]][posMk, 
                  1]
                sebiMk[k, whereisit] <- model.fits[[k]][posMk, 
                  2]
            }
    }
    Ebi <- postprob %*% EbiMk
    Ebimat <- matrix(rep(Ebi, nmod), nrow = nmod, byrow = TRUE)
    SDbi <- sqrt(postprob %*% (sebiMk^2) + postprob %*% ((EbiMk - 
        Ebimat)^2))
    CSDbi <- rep(0, nvar)
    CEbi <- CSDbi
    for (i in (1:ncol(x))) {
        sel <- which[, i]
        if (sum(sel) > 0) {
            cpp <- rbind(postprob[sel]/sum(postprob[sel]))
            CEbi[cox.assign[[i + 1]]] <- as.numeric(cpp %*% EbiMk[sel, 
                cox.assign[[i + 1]]])
            CSDbi[cox.assign[[i + 1]]] <- sqrt(cpp %*% (sebiMk[sel, 
                cox.assign[[i + 1]]]^2) + cpp %*% ((EbiMk[sel, 
                cox.assign[[i + 1]]] - CEbi[cox.assign[[i + 1]]])^2))
        }
    }
    names(output.names) <- var.names
    result <- list(postprob = postprob, label = label, size = size, 
        bic = bic, prior.param = prior.param, prior.model.weights = prior/prior.weight.denom, 
        which = which, probne0 = c(probne0), postmean = as.vector(Ebi), 
        postsd = as.vector(SDbi), condpostmean = CEbi, condpostsd = CSDbi, 
        mle = EbiMk, se = sebiMk, namesx = var.names, reduced = reduced, 
        dropped = dropped, call = cl, n.models = length(postprob), 
        n.vars = length(probne0), nests = length(Ebi), output.names = output.names, 
        assign = cox.assign, factor.type = factor.type)
    class(result) <- "bic.surv"
    result
}
"bic.surv.formula" <-
function (f, data, strict = FALSE, OR = 20, maxCol = 30, prior.param = c(rep(0.5, 
    ncol(x))), OR.fix = 2, nbest = 150, factor.type = TRUE, factor.prior.adjust = FALSE, 
    call = NULL, ...) 
{
    cl <- match.call()
    tms <- terms(f, data = data)
    fmatrix <- attr(tms, "factors")
    tms.order <- attr(tms, "order")
    tms.labels <- attr(tms, "term.labels")
    mm <- model.matrix(tms, data = data)
    assn <- attr(mm, "assign")
    nterms <- max(assn)
    datalist <- eval(attr(tms, "variables"), envir = data)
    nvar <- nrow(fmatrix) - 1
    isvarfac <- rep(NA, times = nvar)
    for (i in 1:nvar) isvarfac[i] <- is.factor(datalist[[i + 
        1]])
    istermfac <- rep(NA, times = nterms)
    for (i in 1:nterms) {
        cterms <- fmatrix[-1, i] == 1
        istermfac[i] <- sum(isvarfac[cterms] == FALSE) == 0
    }
    surv.t.name <- all.vars(f)[1]
    cens.name <- all.vars(f)[2]
    moddata <- data.frame(rep(NA, times = dim(mm)[1]))
    cnames <- NULL
    for (i in 1:nterms) {
        if (istermfac[i]) {
            if (tms.order[i] == 1) {
                moddata <- cbind(moddata, datalist[[i + 1]])
                cnames <- c(cnames, tms.labels[i])
            }
            else {
                sel <- assn == i
                nlev <- sum(sel)
                newfac.index <- (mm[, sel] %*% cbind(1:nlev)) + 
                  1
                facnames <- c("ref", colnames(mm)[sel])
                newfac <- facnames[newfac.index]
                newfac <- factor(newfac)
                moddata <- cbind(moddata, newfac)
                cnames <- c(cnames, paste(tms.labels[i], "..", 
                  sep = ""))
            }
        }
        else {
            sel <- assn == i
            moddata <- cbind(moddata, mm[, sel])
            cnames <- c(cnames, colnames(mm)[sel])
        }
    }
    moddata <- moddata[, -1]
    cnames <- gsub(":", ".", cnames)
    moddata <- cbind(moddata, datalist[[1]][, 1], datalist[[1]][, 
        2])
    colnames(moddata) <- c(cnames, surv.t.name, cens.name)
    nv <- ncol(moddata) - 2
    surv.t <- moddata[, nv + 1]
    cens <- moddata[, nv + 2]
    x <- moddata[, 1:nv]
    bic.surv(x, surv.t, cens, strict = FALSE, OR = OR, maxCol = maxCol, 
        prior.param = prior.param, OR.fix = OR.fix, nbest = nbest, 
        factor.type = factor.type, factor.prior.adjust = factor.prior.adjust, 
        call = cl)
}
"bic.surv.matrix" <-
function (x, surv.t, cens, strict = FALSE, OR = 20, maxCol = 30, 
    prior.param = c(rep(0.5, ncol(x))), OR.fix = 2, nbest = 150, 
    factor.type = TRUE, factor.prior.adjust = FALSE, call = NULL, ...) 
{
    leaps.bs <- function(info, coef, names.arg, nbest = nbest) {
        names.arg <- names.arg
        if (is.null(names.arg)) 
            names.arg <- c(as.character(1:9), LETTERS, letters)[1:ncol(info)]
        if (length(names.arg) < ncol(info)) 
            stop("Too few names")
        bIb <- coef %*% info %*% coef
        kx <- ncol(info)
        maxreg <- nbest * kx
        if (kx < 3) 
            stop("Too few independent variables")
        imeth <- 1
        df <- kx + 1
        Ib <- info %*% coef
        rr <- cbind(info, Ib)
        rr <- rbind(rr, c(Ib, bIb))
        it <- 0
        n.cols <- kx + 1
        nv <- kx + 1
        nf <- 0
        no <- 1e+05
        ib <- 1
        mb <- nbest
        nd <- n.cols
        nc <- 4 * n.cols
        rt <- matrix(rep(0, times = nd * nc), ncol = nc)
        rt[, 1:n.cols] <- rr
        iw <- c(1:(kx + 1), rep(0, times = 4 * nd))
        nw <- length(iw)
        rw <- rep(0, times = 2 * mb * kx + 7 * nd)
        nr <- length(rw)
        t1 <- 2
        s2 <- -1
        ne <- 0
        iv <- 0
        nret <- mb * kx
        Subss <- rep(0, times = nret)
        RSS <- Subss
        ans <- .Fortran("fwleaps", as.integer(nv), as.integer(it), 
            as.integer(kx), as.integer(nf), as.integer(no), as.integer(1), 
            as.double(2), as.integer(mb), as.double(rt), as.integer(nd), 
            as.integer(nc), as.integer(iw), as.integer(nw), as.double(rw), 
            as.integer(nr), as.double(t1), as.double(s2), as.integer(ne), 
            as.integer(iv), as.double(Subss), as.double(RSS), 
            as.integer(nret), PACKAGE = "BMA")
        regid <- ans[[21]]/2
        r2 <- ans[[20]]
        nreg <- sum(regid > 0)
        regid <- regid[1:nreg]
        r2 <- r2[1:nreg]
        which <- matrix(TRUE, nreg, kx)
        z <- regid
        which <- matrix(as.logical((rep.int(z, kx)%/%rep.int(2^((kx - 
            1):0), rep.int(length(z), kx)))%%2), byrow = FALSE, ncol = kx)
        size <- which %*% rep(1, kx)
        label <- character(nreg)
        sep <- if (all(nchar(names.arg) == 1)) 
            ""
        else ","
        for (i in 1:nreg) label[i] <- paste(names.arg[which[i, 
            ]], collapse = sep)
        ans <- list(r2 = r2, size = size, label = label, which = which)
        return(ans)
    }
    factor.names <- function(x) {
        out <- list()
        for (i in 1:ncol(x)) if (is.factor(x[, i])) 
            out[[i]] <- levels(x[, i])
        else out <- c(out, list(NULL))
        attributes(out)$names <- names(x)
        return(out)
    }
    dropcols <- function(x, surv.t, cens, maxCols = 30) {
        vnames <- attributes(x)$names
        nvar <- length(vnames)
        isfac <- rep(FALSE, times = nvar)
        for (i in 1:nvar) isfac[i] <- is.factor(x[, i])
        nlevels <- rep(NA, times = nvar)
        for (i in 1:nvar) if (isfac[i]) 
            nlevels[i] <- length(levels(x[, i]))
        any.dropped <- FALSE
        mm <- model.matrix(terms.formula(~., data = x), data = x)
        designx <- attributes(mm)$assign
        n.designx <- length(designx)
        designx.levels <- rep(1, times = n.designx)
        for (i in 2:n.designx) if (isfac[designx[i]]) 
            designx.levels[i] <- sum(designx[1:i] == designx[i]) + 
                1
        x.coxph <- data.frame(mm[, -1], surv.t = surv.t, cens = cens)
        cox.out <- coxph(Surv(surv.t, cens) ~ ., data = x.coxph, 
            method = "breslow", iter.max = 30)
        while (length(cox.out$coefficients) > maxCol) {
            any.dropped <- TRUE
            dropcox <- drop1(cox.out, test = "Chisq")
            dropped <- which.max(dropcox$"Pr(Chi)"[-1]) + 1
            x.coxph <- x.coxph[, -(dropped - 1)]
            designx.levels <- designx.levels[-dropped]
            designx <- designx[-dropped]
            cox.out <- coxph(Surv(surv.t, cens) ~ ., data = x.coxph, 
                method = "breslow", iter.max = 30)
        }
        remaining.vars <- unique(designx[-1])
        new.nvar <- length(remaining.vars)
        dropped.vars <- vnames[-remaining.vars]
        dropped.levels <- NULL
        ncol.cox <- ncol(x.coxph) - 2
        x.coxph <- x.coxph[-((ncol.cox + 1):(ncol.cox + 2))]
        xx <- data.frame(matrix(rep(NA, times = new.nvar * nrow(x.coxph)), 
            ncol = new.nvar))
        new.names = rep(NA, times = new.nvar)
        for (i in 1:new.nvar) {
            cvar <- remaining.vars[i]
            lvls <- designx.levels[cvar == designx]
            if (isfac[cvar]) {
                if (length(lvls) != length(levels(x[, cvar]))) {
                  newvar <- (as.matrix(x.coxph[, cvar == designx[-1]]) %*% 
                    cbind(lvls - 1)) + 1
                  xx[, i] <- factor(levels(x[, cvar])[newvar])
                  new.names[i] <- vnames[cvar]
                  removed.levels <- levels(x[, cvar])[-c(1, lvls)]
                  dropped.levels <- c(dropped.levels, paste(vnames[cvar], 
                    "_", removed.levels, sep = ""))
                }
                else {
                  xx[, i] <- factor(x[, cvar])
                  new.names[i] <- vnames[cvar]
                }
            }
            else {
                xx[, i] <- x[, cvar]
                new.names[i] <- vnames[cvar]
            }
        }
        dropped <- c(dropped.vars, dropped.levels)
        return(list(mm = xx, any.dropped = any.dropped, dropped = dropped, 
            var.names = new.names, remaining.vars = remaining.vars))
    }
    if (is.null(call)) 
        cl <- match.call()
    else cl <- call
    options(contrasts = c("contr.treatment", "contr.treatment"))
    x <- data.frame(x)
    prior.weight.denom <- 0.5^ncol(x)
    x.omit <- na.omit(x)
    used <- match(row.names(x), row.names(x.omit))
    omitted <- seq(nrow(x))[is.na(used)]
    if (length(omitted) > 0) {
        x <- x.omit
        surv.t <- surv.t[-omitted]
        cens <- cens[-omitted]
        warning(paste("There were ", length(omitted), "records deleted due to NA's"))
    }
    leaps.x <- x
    fac.levels <- rep(1, times = ncol(x))
    fn <- factor.names(x)
    factors <- !all(unlist(lapply(fn, is.null)))
    if (factors) {
        cdf <- cbind.data.frame(y = surv.t, x)
        mm <- model.matrix(formula(cdf), data = cdf)[, -1, drop = FALSE]
        mmm <- data.frame(matrix(mm, nrow = nrow(mm), byrow = FALSE))
        names(mmm) <- dimnames(mm)[[2]]
        output.names <- names(mmm)
        x.coxph <- data.frame(surv.t = surv.t, cens = cens, x)
        cox.out <- coxph(Surv(surv.t, cens) ~ ., data = x.coxph, 
            method = "breslow")
        cox.assign <- cox.out$assign
        for (tempi in length(cox.assign):1) {
            cox.assign[[tempi + 1]] <- cox.assign[[tempi]]
            names(cox.assign)[tempi + 1] <- names(cox.assign)[tempi]
        }
        cox.assign[[1]] <- numeric(0)
        names(cox.assign)[1] <- "(Intercept)"
        fac.levels <- unlist(lapply(cox.assign, length)[-1])
        if (factor.type) {
            for (i in 1:length(names(x))) {
                if (!is.null(fn[[i]])) {
                  nx <- names(x)[i]
                  coefs <- cox.out$coef[cox.out$assign[[i]]]
                  old.vals <- x[, i]
                  new.vals <- c(0, coefs)
                  new.vec <- as.vector(new.vals[match(old.vals, 
                    fn[[i]])])
                  leaps.x[, nx] <- new.vec
                }
            }
        }
        else {
            fac.levels <- rep(1, times = ncol(mmm))
            new.prior <- NULL
            for (i in 1:length(names(x))) {
                addprior <- prior.param[i]
                if (!is.null(fn[[i]])) {
                  k <- length(fn[[i]])
                  if (factor.prior.adjust) 
                    addprior <- rep(1 - (1 - prior.param[i])^(1/(k - 
                      1)), k - 1)
                  else addprior <- rep(prior.param[i], k - 1)
                }
                new.prior <- c(new.prior, addprior)
            }
            prior.param <- new.prior
            x <- leaps.x <- mmm
        }
    }
    xx <- data.frame()
    xx <- dropcols(leaps.x, surv.t, cens, maxCol)
    var.names <- xx$var.names
    remaining <- xx$remaining.vars
    leaps.x <- xx$mm
    reduced <- xx$any.dropped
    dropped <- NULL
    if (reduced) 
        dropped <- xx$dropped
    nvar <- length(x[1, ])
    x <- x[, remaining, drop = FALSE]
    x <- data.frame(x)
    fac.levels <- fac.levels[remaining]
    output.names <- list()
    for (i in 1:length(var.names)) {
        if (is.factor(x[, i])) 
            output.names[[i]] <- levels(x[, i])
        else output.names[[i]] <- NA
    }
    xnames <- names(x)
    names(leaps.x) <- var.names
    x.coxph <- data.frame(surv.t = surv.t, cens = cens, leaps.x)
    cox.out <- coxph(Surv(surv.t, cens) ~ ., data = x.coxph, 
        method = "breslow")
    x.coxph.fac <- data.frame(surv.t = surv.t, cens = cens, x)
    cox.assign <- coxph(Surv(surv.t, cens) ~ ., data = x.coxph.fac, 
        method = "breslow")$assign
    for (tempi in length(cox.assign):1) {
        cox.assign[[tempi + 1]] <- cox.assign[[tempi]]
        names(cox.assign)[tempi + 1] <- names(cox.assign)[tempi]
    }
    cox.assign[[1]] <- numeric(0)
    names(cox.assign)[1] <- "(Intercept)"
    n <- sum(cens)
    if (ncol(leaps.x) >= 3) {
        coef <- cox.out$coef
        info <- solve(cox.out$var)
        a <- leaps.bs(info, coef, names.arg = xnames, nbest = nbest)
        a$r2 <- c(0, a$r2)
        a$size <- c(0, a$size)
        a$label <- c("NULL", a$label)
        a$which <- rbind(rep(FALSE, ncol(x)), a$which)
        nmod <- length(a$size)
        prior.mat <- matrix(rep(prior.param, nmod), nmod, ncol(x), 
            byrow = TRUE)
        prior <- apply(a$which * prior.mat + (!a$which) * (1 - 
            prior.mat), 1, prod)
        bIb <- as.numeric(coef %*% info %*% coef)
        lrt <- bIb - (a$r2 * bIb)
        bic <- lrt + (a$size) * log(n) - 2 * log(prior)
        occam <- bic - min(bic) < 2 * OR.fix * log(OR)
        size <- a$size[occam]
        label <- a$label[occam]
        which <- a$which[occam, , drop = FALSE]
        bic <- bic[occam]
        prior <- prior[occam]
    }
    else {
        nmod <- switch(ncol(x), 2, 4)
        bic <- label <- rep(0, nmod)
        which <- matrix(c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE), nmod, nmod/2)
        size <- c(0, 1, 1, 2)[1:nmod]
        sep <- ","
        prior.mat <- matrix(rep(prior.param, nmod), nmod, ncol(x), 
            byrow = TRUE)
        prior <- apply(which * prior.mat + (!which) * (1 - prior.mat), 
            1, prod)
        for (k in 1:nmod) {
            if (k == 1) 
                label[k] <- "NULL"
            else label[k] <- paste(names(leaps.x)[which[k, ]], 
                collapse = sep)
        }
    }
    model.fits <- as.list(rep(0, length(label)))
    loglik <- rep(0, length(label))
    size <- rep(0, length(label))
    loglik.null <- cox.out$loglik[1]
    for (k in (1:length(label))) {
        if (sum(which[k, ]) != 0) {
            x.coxph <- data.frame(x[, which[k, ], drop = FALSE], 
                surv.t = surv.t, cens = cens)
            cox.out <- coxph(Surv(surv.t, cens) ~ ., data = x.coxph, 
                iter.max = 30, method = "breslow")
            loglik[k] <- cox.out$loglik[2]
            size[k] <- length(cox.out$coef)
            model.fits[[k]] <- matrix(rep(0, 2 * length(cox.out$coef)), 
                ncol = 2)
            model.fits[[k]][, 1] <- cox.out$coef
            model.fits[[k]][, 2] <- sqrt(diag(cox.out$var))
        }
        else {
            loglik[k] <- loglik.null
        }
    }
    bic <- size * log(n) - 2 * (loglik - loglik.null) - 2 * log(prior)
    occam <- bic - min(bic) < 2 * log(OR)
    size <- size[occam]
    label <- label[occam]
    which <- which[occam, , drop = FALSE]
    bic <- bic[occam]
    prior <- prior[occam]
    model.fits <- model.fits[occam]
    postprob <- (exp(-0.5 * (bic - min(bic))))/sum(exp(-0.5 * 
        (bic - min(bic))))
    order.bic <- order(bic, size, label)
    size <- size[order.bic]
    label <- label[order.bic]
    which <- which[order.bic, , drop = FALSE]
    bic <- bic[order.bic]
    prior <- prior[order.bic]
    postprob <- postprob[order.bic]
    model.fits <- model.fits[order.bic]
    nmod <- length(size)
    if (strict & (nmod != 1)) {
        occam <- rep(TRUE, nmod)
        for (k in (2:nmod)) {
            for (j in (1:(k - 1))) {
                which.diff <- which[k, ] - which[j, ]
                if (all(which.diff >= 0)) 
                  occam[k] <- FALSE
            }
        }
        size <- size[occam]
        label <- label[occam]
        which <- which[occam, , drop = FALSE]
        bic <- bic[occam]
        prior <- prior[occam]
        postprob <- postprob[occam]/(sum(postprob[occam]))
        model.fits <- model.fits[occam]
    }
    bic <- bic + 2 * log(prior)
    probne0 <- round(100 * t(which) %*% as.matrix(postprob), 
        1)
    nmod <- length(bic)
    nvar <- max(unlist(cox.assign))
    Ebi <- rep(0, nvar)
    SDbi <- rep(0, nvar)
    EbiMk <- matrix(rep(0, nmod * nvar), nrow = nmod)
    sebiMk <- matrix(rep(0, nmod * nvar), nrow = nmod)
    for (i in (1:ncol(x))) {
        whereisit <- cox.assign[[i + 1]]
        if (any(which[, i])) 
            for (k in (1:nmod)) if (which[k, i] == TRUE) {
                spot <- sum(which[k, (1:i)])
                posMk <- (c(0, cumsum(fac.levels[which[k, ]])) + 
                  1)[spot]
                posMk <- posMk:(posMk + fac.levels[i] - 1)
                EbiMk[k, whereisit] <- model.fits[[k]][posMk, 
                  1]
                sebiMk[k, whereisit] <- model.fits[[k]][posMk, 
                  2]
            }
    }
    Ebi <- postprob %*% EbiMk
    Ebimat <- matrix(rep(Ebi, nmod), nrow = nmod, byrow = TRUE)
    SDbi <- sqrt(postprob %*% (sebiMk^2) + postprob %*% ((EbiMk - 
        Ebimat)^2))
    CSDbi <- rep(0, nvar)
    CEbi <- CSDbi
    for (i in (1:ncol(x))) {
        sel <- which[, i]
        if (sum(sel) > 0) {
            cpp <- rbind(postprob[sel]/sum(postprob[sel]))
            CEbi[cox.assign[[i + 1]]] <- as.numeric(cpp %*% EbiMk[sel, 
                cox.assign[[i + 1]]])
            CSDbi[cox.assign[[i + 1]]] <- sqrt(cpp %*% (sebiMk[sel, 
                cox.assign[[i + 1]]]^2) + cpp %*% ((EbiMk[sel, 
                cox.assign[[i + 1]]] - CEbi[cox.assign[[i + 1]]])^2))
        }
    }
    names(output.names) <- var.names
    result <- list(postprob = postprob, label = label, size = size, 
        bic = bic, prior.param = prior.param, prior.model.weights = prior/prior.weight.denom, 
        which = which, probne0 = c(probne0), postmean = as.vector(Ebi), 
        postsd = as.vector(SDbi), condpostmean = CEbi, condpostsd = CSDbi, 
        mle = EbiMk, se = sebiMk, namesx = var.names, reduced = reduced, 
        dropped = dropped, call = cl, n.models = length(postprob), 
        n.vars = length(probne0), nests = length(Ebi), output.names = output.names, 
        assign = cox.assign, factor.type = factor.type)
    class(result) <- "bic.surv"
    result
}




