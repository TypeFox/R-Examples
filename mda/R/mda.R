mda <-
function(formula = formula(data), data = sys.frame(sys.parent()), 
         subclasses = 3, sub.df = NULL, tot.df = NULL, dimension =
         sum(subclasses) - 1, eps = .Machine$double.eps, iter = 5,
         weights = mda.start(x, g, subclasses, trace, ...), method =
         polyreg, keep.fitted = (n * dimension < 5000), trace = FALSE,
         ...) 
{
    this.call <- match.call()
    m <- match.call(expand.dots = FALSE)
    m[[1]] <- as.name("model.frame")
    m <- m[match(names(m), c("", "formula", "data"), 0)]
    m <- eval(m, sys.frame(sys.parent()))
    Terms <- attr(m, "terms")
    g <- model.extract(m, "response")
    x <- model.matrix(Terms, m)
    if(attr(Terms, "intercept"))x=x[,-1,drop=FALSE]
    dd <- dim(x)
    n <- dd[1]
    m <- dd[2]
    rowmin <- function(mat) {
        ncc <- ncol(mat)
        if (ncc == 1) 
            return(drop(mat))
        rowm <- pmin(mat[, 1], mat[, 2])
        if (ncc == 2) 
            return(rowm)
        else {
            for (j in seq(from = 3, to = ncc)) rowm <- pmin(rowm, 
                mat[, j])
        }
        rowm
    }
    if (length(g) != n) 
        stop("g should have length nrow(x)")
    fg <- factor(g)
    if (inherits(weights, "mda")) {
        if (is.null(weights$weights)) 
            weights <- predict(weights, x, type = "weights", 
                g = fg)
        else weights <- weights$weights
    }
    subclasses <- sapply(weights, ncol)
    prior <- table(fg)
    dim(prior) <- NULL
    prior <- prior/sum(prior)
    cnames <- levels(fg)
    g <- as.numeric(fg)
    J <- length(cnames)
    Assign <- split(seq(sum(subclasses)), rep(seq(J), subclasses))
    names(Assign) <- cnames
    if (!is.null(tot.df)) {
        if (tot.df >= sum(subclasses)) 
            tot.df <- NULL
    }
    if (!is.null(sub.df)) {
        sub.df <- rep(sub.df, length = length(prior))
        sub.df <- pmin(sub.df, subclasses)
        if (all(sub.df == subclasses)) 
            sub.df <- NULL
    }
    for (counter in seq(iter)) {
        fit <- mda.fit(x, g, weights, assign.theta = Assign, 
            Rj = subclasses, sub.df = sub.df, tot.df = tot.df, 
            dimension = dimension, eps = .Machine$double.eps, 
            method = method, trace = trace, ...)
        dmat <- predict.fda(fit, type = "distance")
        sub.prior <- fit$sub.prior
        for (j in seq(J)) {
            TT <- dmat[g == j, Assign[[j]], drop = FALSE]
            TT <- exp(-0.5 * (TT - rowmin(TT)))
            TT <- TT * outer(rep(1, nrow(TT)), sub.prior[[j]])
            weights[[j]][] <- TT/drop(TT %*% rep(1, ncol(TT)))
        }
        pclass <- matrix(1, n, J)
        dmat <- exp(-0.5 * (dmat - rowmin(dmat)))
        for (j in seq(J)) {
            priorj <- sub.prior[[j]]
            ass <- Assign[[j]]
            TT <- dmat[, ass, drop = FALSE] * outer(rep(1, n), priorj)
            TTot <- drop(TT %*% rep(1, length(ass)))
            pclass[, j] <- prior[j] * TTot
        }
        pclass <- pclass/drop(pclass %*% rep(1, J))
        if (trace) 
            cat("Iteration", counter, "\tDeviance(multinomial)", 
                format(round(ll <- llmult(pclass, g), 5)), "\n")
    }
    if (!trace) 
        ll <- llmult(pclass, g)
    if (!keep.fitted) 
        fit$fit$fitted.values <- NULL
    dimnames(pclass) <- list(NULL, names(Assign))
    conf <- confusion(softmax(pclass), fg)
    fit <- c(fit, list(weights = weights, prior = prior, assign.theta = Assign, 
        deviance = ll, confusion = conf, terms = Terms))
    fit$call <- this.call
    fit$sub.df <- sub.df
    fit$tot.df <- tot.df
    class(fit) <- c("mda", "fda")
    fit
}

