
### calculate quad p-values
.MPinv <- function (X, tol = sqrt(.Machine$double.eps)) 
{
    if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
        stop("'X' must be a numeric or complex matrix")
    if (!is.matrix(X)) 
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X)) 
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    Xplus <- if (all(Positive)) 
        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) 
        array(0, dim(X)[2L:1L])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
        t(Xsvd$u[, Positive, drop = FALSE]))
    list(Xplus = Xplus, rank = sum(Positive))
}

.pX2 <- function(lin, exp, cov, pval = TRUE) {
    if (length(lin) == 1) {
        if (cov < .Machine$double.eps) return(c(-Inf, -Inf))
        X2 <- ((lin - exp)^2) / cov
        df <- 1
    } else {
        v <- diag(V <- matrix(cov, ncol = length(lin)))
        ### remove elements with zero variance first
        lin <- as.vector(lin)[v > 0]
        exp <- as.vector(exp)[v > 0]
        V <- V[v > 0, v > 0, drop = FALSE]
        v <- v[v > 0]
        if (length(v) == 0) return(c(-Inf, -Inf))
        tmp <- matrix(lin - exp, ncol = 1)
        Xplus <- .MPinv(V)
        X2 <- crossprod(tmp, Xplus$Xplus) %*% tmp
        df <- Xplus$rank
    }
    
    X2 <- pmax(0, X2) ## Z: X2 may be slightly negative...
    
    if (pval)
        return(c(log(X2), pchisq(X2, df = df, lower.tail = TRUE, 
                                 log.p = TRUE)))
    return(c(log(X2), NA))
}

### calculate max-T p-value
.pmaxT <- function(lin, exp, cov, pval = TRUE) {

    if (length(lin) == 1) {
        if (cov < .Machine$double.eps) return(c(-Inf, -Inf))
        maxT <- abs(lin - exp) / sqrt(cov)
        v <- 1
        V <- matrix(1)
    } else {
        v <- diag(V <- matrix(cov, ncol = length(lin)))
        lin <- as.vector(lin)[v > 0]
        exp <- as.vector(exp)[v > 0]
        V <- V[v > 0, v > 0, drop = FALSE]
        v <- v[v > 0]
        if (length(v) == 0) return(c(-Inf, -Inf))
        maxT <- as.vector(max(abs(lin - exp) / sqrt(v)))
        if (is.na(maxT)) return(c(-Inf, -Inf))
    }
    if (pval) 
        return(c(log(maxT), log(mvtnorm::pmvnorm(lower = rep(-maxT, length(v)),
                                        upper = rep(maxT, length(v)),
                                        sigma = cov2cor(V)))))
    return(c(log(maxT), NA))
}

### surrogate splits
.csurr <- function(split, data, inp, weights, ctrl) {

    ### <FIXME> surrogate splits for multiway splits </FIXME>?
    stopifnot(length(unique(split)) == 2)
    response <- as.factor(split)
    response <- model.matrix(~ response - 1)
    if (ncol(response) == 2) response <- response[, -1, drop = FALSE]
    storage.mode(response) <- "double"

    lin <- .Call("R_LinstatExpCov", data, inp, response, weights)
    ### <FIXME> this is slow, chisq or even teststatistics might be better
    p <- sapply(lin[inp], function(x) do.call(".pX2", x[-1]))
    ### </FIXME>
    colnames(p) <- colnames(data)[inp]
    rownames(p) <- c("statistic", "p.value")
    ### <FIXME> break ties? </FIXME>
    crit <- p["p.value",,drop = TRUE]

    ret <- vector(mode = "list", length = min(sum(inp), ctrl$maxsurrogate))

    for (i in 1L:length(ret)) {
        isel <- which.max(crit)
        isel <- which(inp)[isel]
        x <- data[[isel]]
        sp <- .Call("R_split", x, response, weights, as.integer(0))
        if (any(is.na(sp))) next
        if (length(sp) == 1) {
            ret[[i]] <- partysplit(as.integer(isel), breaks = sp, index = 1L:2L)
        } else {
            ret[[i]] <- partysplit(as.integer(isel), index = sp)
        }
        tmp <- kidids_split(ret[[i]], data, obs = weights > 0)
        tmps <- split[weights > 0]
        tab <- table(tmp, tmps)
        if (tab[1, 1] < tab[1, 2]) {
            indx <- ret[[i]]$index
            ret[[i]]$index[indx == 1] <- 2
            ret[[i]]$index[indx == 2] <- 1
        }
        crit[which.max(crit)] <- -Inf
    }
    ret <- ret[!sapply(ret, is.null)]
    if (length(ret) == 0) ret <- NULL
    return(ret)
}

### set up new node for conditional inference tree
.cnode <- function(id = 1, data, response, inputs, weights, ctrl, cenv = NULL) {

    if (is.null(cenv)) {
        cenv <- new.env()
        depth <- 0
    } else {
        depth <- get("depth", envir = cenv)
        if (depth >= ctrl$maxdepth)
            return(partynode(as.integer(id)))
    }
    weights <- as.integer(weights)
    if (sum(weights) < ctrl$minsplit) return(partynode(as.integer(id)))
    if (id > 1 && ctrl$stump) return(partynode(as.integer(id)))

    inp <- inputs
    if (ctrl$mtry < Inf) {
        mtry <- min(sum(inp), ctrl$mtry)
        ### sum(inp) == 1 will lead to sample.int instead of sample; see example(sample)
        resample <- function(x, ...) x[sample.int(length(x), ...)]
        s <- resample(which(inp), mtry)
        inp <- logical(length(inp))
        inp[s] <- TRUE
    } 

    response_arg <- response
    if (is.function(response))
        response <- response(data, weights)

    lin <- .Call("R_LinstatExpCov", data, inp, response, weights)
    ### (potentially) parallel computation of criterion
    p <- simplify2array(ctrl$applyfun(lin[inp], function(x) 
        do.call(ctrl$cfun, x[-1])))
    crit <- p[1,,drop = TRUE]
    ### crit is maximised, but there might be ties
    ties <- which(abs(crit - max(crit)) < .Machine$double.eps)
    if (length(ties) > 1) {
        ### add a small value (< 1/1000) to crit derived from order of teststats
        crit[ties] <- crit[ties] + order(p["statistic", ties]) / (sum(ties) * 1000)
    }
    p <- p[-1,,drop = FALSE]
    colnames(p) <- colnames(data)[inp]

    mb <- ctrl$minbucket
    mp <- ctrl$minprob
    storage.mode(mb) <- "integer"

    ### format p values
    fmP <- function(p) {
        if (all(is.na(p["p.value",]))) return(NA)
        1 - exp(p["p.value",])
    }

    count <- 1
    thissplit <- NULL
    while(count <= ctrl$splittry) {
        if (any(crit > ctrl$mincriterion)) {
            isel <- iselp <- which.max(crit)
            isel <- which(inp)[isel]
        } else {
            return(partynode(as.integer(id), 
                             info = list(criterion = p,
                                         p.value = min(fmP(p), na.rm = TRUE))))
        }
        x <- data[[isel]]
        swp <- ceiling(sum(weights) * mp)
        if (mb < swp) mb <- as.integer(swp)

        if ((ctrl$multiway && ctrl$maxsurrogate == 0) && is.factor(x)) {
            if (all(table(x[rep(1L:length(x), weights)]) > mb)) {
                thissplit <- partysplit(as.integer(isel), index = 1L:nlevels(x))
                break()
            }
        } else {
            sp <- .Call("R_split", x, response, weights, mb)
            if (!any(is.na(sp))) {
                if (length(sp) == 1) {
                    thissplit <- partysplit(as.integer(isel), breaks = sp)
                } else {
                    ### deal with empty levels -> NA in sp
                    if (is.factor(x)) 
                        sp[table(rep(x, weights)) == 0] <- NA
                    thissplit <- partysplit(as.integer(isel), index = sp)
                }
                break()
            }
        }
        crit[which.max(crit)] <- -Inf
        count <- count + 1
    }
    if (is.null(thissplit))
        return(partynode(as.integer(id), 
                         info = list(criterion = p,
                                     p.value = min(fmP(p), na.rm = TRUE))))           

    ret <- partynode(as.integer(id))
    ret$split <- thissplit
    ret$info <- list(criterion = p, p.value = fmP(p)[iselp])
    thissurr <- NULL
    kidids <- kidids_node(ret, data)
    prob <- prop.table(table(kidids))
    if (ctrl$majority)  ### go with majority
        prob <- numeric(0) + 1L:length(prob) %in% which.max(prob)
    ret$prob <- prob

    if (ctrl$maxsurrogate > 0) {
        inp <- inputs
        inp[isel] <- FALSE
        w <- weights
        xna <- is.na(x)
        w[xna] <- 0L
        ret$surrogates <- .csurr(kidids, data, inp, w, ctrl)
        kidids[xna] <- kidids_node(ret, data, obs = xna)
    }

    kids <- vector(mode = "list", length = max(kidids)) ## Z: was 1:max(kidids)
    nextid <- id + 1
    for (k in 1L:max(kidids)) {
        w <- weights
        w[kidids != k] <- 0
        assign("depth", depth + 1, envir = cenv)
        kids[[k]] <- .cnode(nextid, data, response_arg, inputs, w, ctrl, cenv)
        nextid <- max(nodeids(kids[[k]])) + 1
    }
    ret$kids <- kids

    return(ret)
}

ctree_control <- function(teststat = c("quad", "max"),
    testtype = c("Bonferroni", "Univariate", "Teststatistic"),
    mincriterion = 0.95, minsplit = 20L, minbucket = 7L, minprob = 0.01,
    stump = FALSE, maxsurrogate = 0L, mtry = Inf, maxdepth = Inf, 
    multiway = FALSE, splittry = 2L, majority = FALSE,
    applyfun = NULL, cores = NULL) {

    teststat <- match.arg(teststat)
    testtype <- match.arg(testtype)

    ## apply infrastructure for determining split points
    if (is.null(applyfun)) {
        applyfun <- if(is.null(cores)) {
            lapply
        } else {
            function(X, FUN, ...) 
                parallel::mclapply(X, FUN, ..., mc.cores = cores)
        }
    }


    list(teststat = teststat,
         testtype = testtype, mincriterion = log(mincriterion),
         minsplit = minsplit, minbucket = minbucket, 
         minprob = minprob, stump = stump, mtry = mtry, 
         maxdepth = maxdepth, multiway = multiway, splittry = splittry,
         maxsurrogate = maxsurrogate, majority = majority, 
         applyfun = applyfun)
}

ctree <- function(formula, data, weights, subset, na.action = na.pass, 
                  control = ctree_control(...), ytrafo = NULL, 
                  scores = NULL, ...) {

    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    
    ### only necessary for extended model formulae 
    ### e.g. multivariate responses
    formula <- Formula::Formula(formula)
    mf$formula <- formula
    mf$drop.unused.levels <- FALSE
    mf$na.action <- na.action
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    response <- names(Formula::model.part(formula, mf, lhs = 1))
    weights <- model.weights(mf)
    dat <- mf[, colnames(mf) != "(weights)"]
    if (!is.null(scores)) {
        for (n in names(scores)) {
            sc <- scores[[n]]
            if (is.ordered(dat[[n]]) && 
                nlevels(dat[[n]]) == length(sc)) {
                attr(dat[[n]], "scores") <- as.numeric(sc)
            } else {
                warning("scores for variable ", sQuote(n), " ignored")
            }
        }
    }

    if (is.null(weights))
        weights <- rep(1, nrow(mf))
    storage.mode(weights) <- "integer"

    nvar <- sum(!(colnames(dat) %in% response))

    control$cfun <- function(...) {
        if (control$teststat == "quad")
            p <- .pX2(..., pval = (control$testtype != "Teststatistic"))
        if (control$teststat == "max")
            p <- .pmaxT(..., pval = (control$testtype != "Teststatistic"))
        names(p) <- c("statistic", "p.value")

        if (control$testtype == "Bonferroni")
            p["p.value"] <- p["p.value"] * min(nvar, control$mtry)
        crit <-  p["statistic"]
        if (control$testtype != "Teststatistic")
        crit <- p["p.value"]
        c(crit, p)
    }

    tree <- .ctree_fit(dat, response, weights = weights, ctrl = control, 
                       ytrafo = ytrafo)

    fitted <- data.frame("(fitted)" = fitted_node(tree, dat), 
                         "(weights)" = weights,
                         check.names = FALSE)
    fitted[[3]] <- dat[, response, drop = length(response) == 1]
    names(fitted)[3] <- "(response)"
    ret <- party(tree, data = dat, fitted = fitted, 
                 info = list(call = match.call(), control = control))
    class(ret) <- c("constparty", class(ret))

    ### doesn't work for Surv objects
    # ret$terms <- terms(formula, data = mf)
    ret$terms <- terms(mf)
    ### need to adjust print and plot methods
    ### for multivariate responses
    ### if (length(response) > 1) class(ret) <- "party"
    return(ret)
}

.ctree_fit <- function(data, response, weights = NULL,
                       ctrl, ytrafo = NULL) {

    inputs <- !(colnames(data) %in% response)

    if (is.null(weights))
        weights <- rep(1, nrow(data))
    storage.mode(weights) <- "integer"

    ### <FIXME> this interface has to change; we need to be
    ### closer to mob() with y ~ x | z formulae and probably a `fit'
    ### argument </FIXME>
    if (!is.function(ytrafo)) {
        infl <- .y2infl(data, response, ytrafo = ytrafo)
    } else {
        infl <- ytrafo ### will be updated with weights in every node
    }

    tree <- .cnode(1L, data, infl, inputs, weights, ctrl)

    return(tree)
}

#.logrank_trafo <- function(x, ties.method = c("logrank", "HL")) {
#    ties.method <- match.arg(ties.method)
#    time <- x[,1]
#    event <- x[,2]
#    n <- length(time)
#    ot <- order(time, event)
#    rt <- rank(time, ties.method = "max")
#    mt <- rank(time, ties.method = "min") - 1
#    fact <- switch(ties.method, "logrank" = event / (n - mt),
#                                "HL" = event/(n - rt + 1)
#                  )   
#    event - cumsum(fact[ot])[rt]
#}

.logrank_trafo <- function(...)
    return(coin::logrank_trafo(...))

### convert response y to influence function h(y)
.y2infl <- function(data, response, ytrafo = NULL) {

    if (length(response) == 1) {
        if (!is.null(ytrafo[[response]])) {
            yfun <- ytrafo[[response]]
            rtype <- "user-defined"
        } else {
            rtype <- class(data[[response]])[1]
            if (rtype == "integer") rtype <- "numeric"
        }
        response <- data[[response]]

        infl <- switch(rtype,
            "user-defined" = yfun(response),
            "factor" = { 
                X <- model.matrix(~ response - 1)
                if (nlevels(response) > 2) return(X)
                return(X[,-1, drop = FALSE])
            },
            "ordered" = {
                sc <- attr(response, "scores")
                if (is.null(sc)) sc <- 1L:nlevels(response)
                sc <- as.numeric(sc)
                return(sc[as.integer(response)])
            },
            "numeric" = response,
            "Surv" = .logrank_trafo(response)
        )
    } else {
        ### multivariate response
        infl <- lapply(response, .y2infl, data = data)
        tmp <- do.call("cbind", infl)
        attr(tmp, "assign") <- rep(1L:length(infl), sapply(infl, NCOL))
        infl <- tmp
    }
    storage.mode(infl) <- "double"
    return(infl)
}

sctest.constparty <- function(object, node = NULL, ...)
{

    if(is.null(node)) {
        ids <- nodeids(object, terminal = FALSE) ### all nodes
    } else {
        ids <- node
    }

    rval <- nodeapply(object, ids, function(n) {
        crit <- info_node(n)$criterion
        if (is.null(crit)) return(NULL)
        ret <- exp(crit)
        ret["p.value",] <- 1 - ret["p.value",]
        ret
    })
    names(rval) <- ids
    if(length(ids) == 1L)
        return(rval[[1L]])
    return(rval)
}
