bicreg <-
function (x, y, wt = rep(1, length(y)), strict = FALSE, OR = 20, 
    maxCol = 31, drop.factor.levels = TRUE, nbest = 150) 
{
    dropcols <- function(x, y, wt, maxCols = 31) {
        x1.ldf <- data.frame(x, y = y)
        temp.wt <- wt
        lm.out <- lm(y ~ ., data = x1.ldf, weights = temp.wt)
        form.vars <- all.vars(formula(lm.out))[-1]
        any.dropped <- FALSE
        dropped.which <- NULL
        while (length(lm.out$coefficients) > maxCol) {
            any.dropped <- TRUE
            droplm <- drop1(lm.out, test = "none")
            dropped <- row.names(droplm)[which.min(droplm$RSS[-1]) + 
                1]
            dropped.index <- match(dropped, form.vars)
            form.vars <- form.vars[-dropped.index]
            formla <- formula(paste("y", "~", paste(form.vars, 
                collapse = " + "), sep = " "))
            lm.out <- lm(formla, data = x1.ldf, weights = temp.wt)
            dropped.which <- c(dropped.which, dropped)
        }
        new.var.names <- names(lm.out$coefficients)
        return(list(mm = model.matrix(lm.out)[, -1, drop = FALSE], 
            any.dropped = any.dropped, dropped = dropped.which, 
            var.names = new.var.names))
    }
    cl <- match.call()
    x <- data.frame(x)
    if (is.null(dimnames(x))) 
        dimnames(x) <- list(NULL, paste("X", 1:ncol(x), sep = ""))
    y <- as.numeric(y)
    options(contrasts = c("contr.treatment", "contr.treatment"))
    xnames <- input.names <- dimnames(x)[[2]]
    x2 <- na.omit(data.frame(x))
    used <- match(row.names(data.frame(x)), row.names(x2))
    omitted <- seq(nrow(x))[is.na(used)]
    if (length(omitted) > 0) {
        wt <- wt[-omitted]
        x <- x2
        y <- y[-omitted]
        warning(paste("There were ", length(omitted), "records deleted due to NA's"))
    }
    if (drop.factor.levels) {
        cdf <- cbind.data.frame(y = y, x)
        mm <- model.matrix(formula(cdf), data = cdf)[, -1, drop = FALSE]
        x <- mm
    }
    xx <- dropcols(x, y, wt, maxCol)
    xnames <- xx$var.names[-1]
    x <- xx$mm
    reduced <- xx$any.dropped
    dropped <- NULL
    if (reduced) 
        dropped <- xx$dropped
    nvar <- length(x[1, ])
    if (nvar > 2) {
#       a <- leaps(x, y, wt = wt, method = "r2", names = dimnames(x)[[2]], 
#            strictly.compatible = FALSE, nbest = nbest)
        a <- regsubsets(x, y, weights = wt, nbest = nbest,
                        nvmax = ncol(x), method = "exhaustive",
                        really.big = TRUE)
        a    <- summary(a)
        size <- apply(a$which,1,sum)
        names(size) <- NULL
        size <- c(1,size)
        a$which <- a$which[,-1, drop=FALSE]
        a$r2 <- a$rsq
        a$r2 <- pmin(pmax(0, a$r2), 0.999)
        x.lm <- cbind.data.frame(y = y, as.data.frame(x[, a$which[2, 
            , drop = FALSE]]), w = wt)
        lm.fix <- lm(y ~ . - w, weights = wt, data = x.lm)
        r2.fix <- summary(lm.fix)$r.sq
        N <- ncol(x)
        magic <- N * log(1 - a$r2[2]) - N * log(1 - r2.fix)
        a$r2 <- 1 - (1 - a$r2) * exp(-magic/N)
        r2 <- round(c(0, a$r2) * 100, 3)
        which <- rbind(rep(FALSE, ncol(x)), a$which)
        templabs <- t(matrix(rep(colnames(which), times = nrow(which)), 
            ncol = nrow(which)))
        templabs[!which] <- ""
        label <- apply(templabs, 1, paste, collapse = "")
        label[1] <- "NULL"
    }
    else {
        r2 <- bic <- NULL
        nmod <- switch(ncol(x), 2, 4)
        bic <- label <- rep(0, nmod)
        model.fits <- as.list(rep(0, nmod))
        which <- matrix(c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, 
            TRUE, TRUE), nmod, nmod/2)
        size <- c(1, 2, 2, 3)[1:nmod]
        sep <- if (all(nchar(dimnames(x)[[2]]) == 1)) 
            ""
        else ","
        for (k in 1:nmod) {
            if (k == 1) {
                label[k] <- "NULL"
                lm1 <- lm(y ~ 1, weights = wt)
            }
            else {
                label[k] <- paste(dimnames(x)[[2]][which[k, ]], 
                  collapse = sep)
                x.lm <- cbind.data.frame(y = y, x = x[, which[k, 
                  , drop = FALSE]], wt = wt)
                lm1 <- lm(y ~ . - wt, data = x.lm, weights = wt)
            }
            r2[k] <- summary(lm1)$r.sq * 100
        }
    }
    n <- length(y)
    if (any((1 - r2/100) <= 0)) 
        stop("a model is perfectly correlated with the response")
    bic <- n * log(1 - r2/100) + (size - 1) * log(n)
    occam <- bic - min(bic) < 2 * log(OR)
    r2 <- r2[occam]
    size <- size[occam]
    label <- label[occam]
    which <- which[occam, , drop = FALSE]
    bic <- bic[occam]
    mbic <- bic - max(bic)
    postprob <- exp(-0.5 * mbic)/sum(exp(-0.5 * mbic))
    postprob[is.na(postprob)] <- 1
    order.bic <- order(bic, size, label)
    r2 <- r2[order.bic]
    size <- size[order.bic]
    label <- label[order.bic]
    which <- which[order.bic, , drop = FALSE]
    bic <- bic[order.bic]
    postprob <- postprob[order.bic]
    if (strict) {
        nmod <- length(bic)
        if (nmod > 1) {
            occam <- rep(TRUE, nmod)
            for (k in (2:nmod)) {
                for (j in (1:(k - 1))) {
                  which.diff <- which[k, ] - which[j, ]
                  if (all(which.diff >= 0)) 
                    occam[k] <- FALSE
                }
            }
            r2 <- r2[occam]
            size <- size[occam]
            label <- label[occam]
            nmod <- sum(occam)
            which <- which[occam, , drop = FALSE]
            bic <- bic[occam]
            postprob <- postprob[occam]
            postprob <- postprob/sum(postprob)
        }
    }
    probne0 <- round(100 * t(which) %*% as.matrix(postprob), 
        1)
    nmod <- length(bic)
    model.fits <- as.list(rep(0, nmod))
    for (k in (1:nmod)) {
        if (sum(which[k, ]) != 0) {
            model.fits[[k]] <- ls.print(lsfit(x[, which[k, ], 
                drop = FALSE], y, wt = wt), print.it = FALSE)$coef.table[[1]]
        }
        else model.fits[[k]] <- ls.print(lsfit(rep(1, length(y)), 
            y, wt = wt, intercept = FALSE), print.it = FALSE)$coef.table[[1]]
    }
    Ebi <- rep(0, (nvar + 1))
    SDbi <- rep(0, (nvar + 1))
    CEbi <- Ebi
    CSDbi <- SDbi
    EbiMk <- matrix(rep(0, nmod * (nvar + 1)), nrow = nmod)
    sebiMk <- matrix(rep(0, nmod * (nvar + 1)), nrow = nmod)
    for (i in 1:(nvar + 1)) {
        if ((i == 1) || (sum(which[, (i - 1)] != 0))) {
            for (k in (1:nmod)) {
                if ((i == 1) || (which[k, (i - 1)] == TRUE)) {
                  if (i == 1) 
                    pos <- 1
                  else pos <- 1 + sum(which[k, (1:(i - 1))])
                  EbiMk[k, i] <- model.fits[[k]][pos, 1]
                  sebiMk[k, i] <- model.fits[[k]][pos, 2]
                }
            }
            Ebi[i] <- as.numeric(sum(postprob * EbiMk[, i]))
            SDbi[i] <- sqrt(postprob %*% (sebiMk[, i]^2) + postprob %*% 
                ((EbiMk[, i] - Ebi[i])^2))
            if (i == 1) {
                CEbi[i] <- Ebi[i]
                CSDbi[i] <- SDbi[i]
            }
            else {
                sel <- which[, i - 1]
                cpp <- postprob[sel]/sum(postprob[sel])
                CEbi[i] <- as.numeric(sum(cpp * EbiMk[sel, i]))
                CSDbi[i] <- sqrt(cpp %*% (sebiMk[sel, i]^2) + 
                  cpp %*% ((EbiMk[sel, i] - CEbi[i])^2))
            }
        }
    }
    dimnames(which) <- list(NULL, colnames(x))
    dimnames(EbiMk) <- dimnames(sebiMk) <- list(NULL, c("(Intercept)", 
        colnames(x)))
    postmean <- apply(postprob * EbiMk, 2, sum)
    fittedValues <- function(coef, x) cbind(1, x) %*% coef
    residualVariance <- function(coef, x, y) sum((y - fittedValues(coef, 
        x))^2)/(length(y) - length(coef))
    resvar <- apply(EbiMk, 1, residualVariance, x = x, y = y)
    result <- list(postprob = postprob, namesx = xnames, label = label, 
        r2 = r2, bic = bic, size = (size - 1), which = which, 
        probne0 = c(probne0), postmean = postmean, residvar = resvar, 
        postsd = SDbi, condpostmean = CEbi, condpostsd = CSDbi, 
        ols = EbiMk, mle = EbiMk, se = sebiMk, reduced = reduced, 
        dropped = dropped, input.names = input.names, call = cl, 
        n.models = length(postprob), n.vars = length(probne0))
    class(result) <- "bicreg"
    result
}
