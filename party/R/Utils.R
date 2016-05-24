
# $Id: Utils.R 537 2014-06-26 12:09:22Z thothorn $

### Wrapper for functions defined in ./src/Utilsc

qsvd <- function(x) {
    if (!is.matrix(x) || ncol(x) != nrow(x))
        stop(sQuote("x"), " is not a quadratic matrix")

    svdmem <- new("svd_mem", ncol(x)) 
    dummy <- .Call("R_svd", x, svdmem, PACKAGE = "party")
    rm(dummy)
    return(list(u = svdmem@u, vt = svdmem@v, d = svdmem@s))
}

MPinv <- function(x, tol = sqrt(.Machine$double.eps)) {
    if (!is.matrix(x) || ncol(x) != nrow(x))
        stop(sQuote("x"), " is not a quadratic matrix")

    svdmem <- new("svd_mem", ncol(x))
    RET <- .Call("R_MPinv", x, tol, svdmem, PACKAGE = "party")
    return(RET@MPinv)
}

### wrapper for functions defined in ./src/Splits.c

Split <- function(x, y, weights, splitctrl) {

    if (is.factor(y))
        ym <- sapply(levels(y), function(l) as.numeric(y == l))
     else 
        ym <- matrix(y, ncol = 1)
    storage.mode(ym) <- "double"

    if (is.factor(x)) {
        xm <- sapply(levels(x), function(l) as.numeric(x == l))
        storage.mode(xm) <- "double"
        xc <- as.numeric(x)
        storage.mode(xc) <- "integer"
        lecxy <- new("LinStatExpectCovar", ncol(xm), ncol(ym))
        lec <- new("LinStatExpectCovar", as.integer(1), ncol(ym))
        eci <- ExpectCovarInfluence(ym, weights)
        split <- .Call("R_splitcategorical", xm, xc, ym, weights, lec, lecxy,
                       eci, splitctrl, PACKAGE = "party")
    } else {
        ox <- order(x)
        storage.mode(ox) <- "integer"
        xm <- matrix(x, ncol = 1)
        storage.mode(xm) <- "double"
        lec <- new("LinStatExpectCovar", as.integer(1), ncol(ym))
        eci <- ExpectCovarInfluence(ym, weights)
        split <- .Call("R_split", xm, ym, weights, ox, lec,
                       eci, splitctrl, PACKAGE = "party")
    }
    split
}

### Wrapper for functions defined in ./src/TestStatistic.c

maxabsTestStatistic <- function(t, mu, Sigma, tol = sqrt(.Machine$double.eps)) {
    storage.mode(t) <- "double"
    storage.mode(mu) <- "double"
    storage.mode(Sigma) <- "double"
    storage.mode(tol) <- "double"
    
    if (length(t) != length(mu) || length(t) != nrow(Sigma)) 
        stop("dimensions don't match")
        
    .Call("R_maxabsTestStatistic", t, mu, Sigma, tol, PACKAGE = "party")
}

quadformTestStatistic <- function(t, mu, Sigma, 
    tol = sqrt(.Machine$double.eps)) {

    storage.mode(t) <- "double"
    storage.mode(mu) <- "double"
    storage.mode(Sigma) <- "double"
    storage.mode(tol) <- "double"
    
    if (length(t) != length(mu) || length(t) != nrow(Sigma)) 
        stop("dimensions don't match")
        
    SigmaPlus <- MPinv(Sigma, tol = tol)
    .Call("R_quadformTestStatistic", t, mu, SigmaPlus, PACKAGE = "party")
}


### Wrapper for functions defined in ./src/LinearStatistic.c

LinearStatistic <- function(x, y, weights) {
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    .Call("R_LinearStatistic", x, y, weights, PACKAGE = "party")
}

ExpectCovarInfluence <- function(y, weights) {
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    .Call("R_ExpectCovarInfluence", y, weights, PACKAGE = "party")
}

ExpectCovarLinearStatistic <- function(x, y, weights) {
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    expcovinf <- ExpectCovarInfluence(y, weights)
    .Call("R_ExpectCovarLinearStatistic", x, y, weights, expcovinf,
          PACKAGE = "party")
}

PermutedLinearStatistic <- function(x, y, indx, perm) {
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"

    if (any(indx < 1 || indx > nrow(y)))
        stop("wrong indices")

    if (any(perm < 1 || perm > nrow(y)))
        stop("wrong indices")

    # C indexing 
    indx <- indx - 1
    perm <- perm - 1
    storage.mode(indx) <- "integer"
    storage.mode(perm) <- "integer"
    .Call("R_PermutedLinearStatistic", x, y, indx, perm, 
          PACKAGE = "party")
}

### Median Survival Time, see survival:::print.survfit

mst <- function(x) {
    minmin <- function(y, xx) {

        ### don't complain about Inf 
        ww <- getOption("warn")
        on.exit(options(warn = ww))
        options(warn = -1)

        if (any(!is.na(y) & y==.5)) {
            if (any(!is.na(y) & y <.5))
                .5*(min(xx[!is.na(y) & y==.5]) + min(xx[!is.na(y) & y<.5]))
            else
                .5*(min(xx[!is.na(y) & y==.5]) + max(xx[!is.na(y) & y==.5]))
        } else   min(xx[!is.na(y) & y<=.5])
    }
    med <- minmin(x$surv, x$time)
    return(med)
}

dostep <- function(x, y) {

    ### create a step function based on x, y coordinates
    ### modified from `survival:print.survfit'
    if (is.na(x[1] + y[1])) {
        x <- x[-1]
        y <- y[-1]
    }
    n <- length(x)
    if (n > 2) {  
        # replace verbose horizonal sequences like
        # (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
        # with (1, .2), (3, .1).  They are slow, and can smear the looks
        # of the line type.
        dupy <- c(TRUE, diff(y[-n]) !=0, TRUE)
        n2 <- sum(dupy)

        #create a step function
        xrep <- rep(x[dupy], c(1, rep(2, n2-1)))
        yrep <- rep(y[dupy], c(rep(2, n2-1), 1))
        RET <- list(x = xrep, y = yrep)
    } else {
        if (n == 1) {
            RET <- list(x = x, y = y)
        } else {
            RET <- list(x = x[c(1,2,2)], y = y[c(1,1,2)])
        }
    }
    return(RET)
}

### Normalized mutual information, Stehl & Gosh (2002) JMLR
nmi <- function(x, y)
{
    x <- table(x, y)
    x <- x / sum(x)                     # ???

    m_x <- rowSums(x)
    m_y <- colSums(x)
    y <- outer(m_x, m_y)

    i <- which((x > 0) & (y > 0))
    out <- sum(x[i] * log(x[i] / y[i]))
    e_x <- sum(m_x * log(ifelse(m_x > 0, m_x, 1)))
    e_y <- sum(m_y * log(ifelse(m_y > 0, m_y, 1)))

    out / sqrt(e_x * e_y)
}

### check if two objects are identical and print differences else
isequal <- function(a, b) {
    attributes(a) <- NULL
    attributes(b) <- NULL
    if (!isTRUE(all.equal(a, b))) {
        print(a, digits = 10)
        print(b, digits = 10)
        return(FALSE)
    } else {
        return(TRUE)
    }
}

mysurvfit <- function(y, weights, ...) {

    stopifnot(extends(class(y), "Surv"))
    ### see comment on weights and subset in ?survfit
    y <- y[weights > 0,]
    weights <- weights[weights > 0]
    return(survfit(y ~ 1, weights = weights, ...))
}

R_get_nodeID <- function(tree, inputs, mincriterion)
    .Call("R_get_nodeID", tree, inputs, 0.0, -1L, PACKAGE = "party")

R_getpredictions <- function(tree, where)
    .Call("R_getpredictions", tree, where, PACKAGE = "party")

R_remove_weights <- function(tree, remove)
    .Call("R_remove_weights", tree, remove, package = "party")

R_modify_response <- function(y, responses)
    .Call("R_modify_response", as.double(y), responses,
          PACKAGE = "party")

R_TreeGrow <- function(object, weights, fitmem, ctrl, where)
    .Call("R_TreeGrow", object, weights, fitmem, ctrl,
          where, PACKAGE = "party")

copyslots <- function(source, target) {
    slots <- names(getSlots(class(source)))
    slots <- slots[(slots %in% names(getSlots(class(target))))]
    if (length(slots) == 0) 
        stop("no common slots to copy to")
    for (s in slots)
        eval(parse(text = paste("target@", s, " <- source@", s)))
    return(target)
}
