.capitalize <- function(x) {
    return(paste(toupper(substring(x, 1, 1)), substring(x, 2), sep=''))
}


.hedge <- function(x, hedgeVals) {
    stopifnot(dim(x) == 1)
    stopifnot(length(hedgeVals) == 3)

    if (length(x) <= 0) {
        return(numeric(0))
    }
    if (length(x) > 1) {
        return(sapply(x, .hedge, hedgeVals))
    }

    if (x <= hedgeVals[1]) {
        return(0)
    }
    if (x <= hedgeVals[2]) {
        t <- x - hedgeVals[1]
        t <- t * t / ((hedgeVals[2] - hedgeVals[1]) * (hedgeVals[3] - hedgeVals[1]))
        return(t)
    }
    if (x <= hedgeVals[3]) {
        t <- hedgeVals[3] - x
        t <- 1 - t * t / ((hedgeVals[3] - hedgeVals[2]) * (hedgeVals[3] - hedgeVals[1]))
        return(t)
    }
    return(1)
}


.hedgize3 <- function(x, 
                      context,
                      name,
                      type=c('sm', 'me', 'bi'),
                      allowed=c("ex", "si", "ve", "ml", "ro", "qr", "vr")) {
    type <- match.arg(type)
    hedges <- .fset3[.fset3[, type], ]
    allowed <- match.arg(allowed, several.ok=TRUE)
    allowed <- c("", allowed)
    allowed <- intersect(allowed, hedges$short)


    low <- context[1]
    center <- context[2]
    big <- context[3]

    if (length(x) <= 0) {
        horizon <- numeric(0)
    } else if (type == 'sm') {
        horizon <- ifelse(low <= x & x < center,
                          (center - x) / (center - low),
                          0)
    } else if (type == 'me') {
        horizon <- ifelse(low < x & x <= center,
                          (x - low) / (center - low),
                          ifelse(center <= x & x < big,
                                 (x - big) / (center - big),
                                 0))
    } else if (type == 'bi') {
        horizon <- ifelse(center < x & x <= big,
                          (x - center) / (big - center),
                          0)
    } else {
        horizon <- NaN
    }

    res <- NULL
    if (length(horizon) <= 0) {
        res <- matrix(0, nrow=0, ncol=length(allowed))
    } else {
        res <- sapply(allowed, function(sh) {
            .hedge(horizon, as.vector(as.matrix((hedges[hedges$short==sh, c('p1', 'p2', 'p3')]))))
        })
    }
    if (!is.matrix(res)) {
        res <- t(as.matrix(res))
    }

    colnames(res) <- paste(.capitalize(allowed), .capitalize(type), '.', name, sep='')
    rownames(res) <- NULL

    vars <- rep(name, ncol(res))
    names(vars) <- colnames(res)

    specs <- matrix(0, nrow=nrow(.fset3), ncol=nrow(.fset3))
    specs[row(specs) > col(specs)] <- 1
    colnames(specs) <- paste(.capitalize(.fset3$short), .capitalize(type), '.', name, sep='')
    rownames(specs) <- colnames(specs)
    specs <- specs[colnames(res), colnames(res), drop=FALSE]

    return(fsets(res, vars=vars, specs=specs))
}


lcut3.numeric <- function(x, 
                          context=NULL,
                          defaultCenter=0.5,
                          atomic=c("sm", "me", "bi"),
                          hedges=c("ex", "si", "ve", "ml", "ro", "qr", "vr"),
                          name=NULL,
                          parallel=FALSE,
                          ...) {
    if (!is.vector(x) || !is.numeric(x)) {
        stop("'x' is not a numeric vector")
    }
    if (is.null(context)) {
        lo <- min(x)
        hi <- max(x)
        context <-c(lo, (hi-lo) * defaultCenter + lo, hi)
    }
    if (length(context) != 3 || context[1] != min(context) || context[3] != max(context)) {
        stop("'context' must be vector with 3 values (lo, med, hi) where lo <= med <= hi")
    }
    #if (context[1] >= context[3]) {
        #stop("'context[1]' must be lower than 'context[3]'")
    #}
    if (is.null(name)) {
        stop("If 'x' is numeric vector then 'name' must not be NULL")
    }
    if (is.null(defaultCenter) || defaultCenter < 0 || defaultCenter > 1) {
        stop("'defaultCenter' must be a number in the interval [0, 1]")
    }
    if (!is.logical(parallel) || length(parallel) != 1) {
        stop("'parallel' must be either TRUE or FALSE")
    }

    x[x < context[1]] <- context[1]
    x[x > context[3]] <- context[3]

    hedges <- match.arg(hedges, several.ok=TRUE)
    atomic <- match.arg(atomic, several.ok=TRUE)

    if (length(atomic) <= 0) {
        stop("'atomic' must not be empty")
    }

    sm <- NULL
    me <- NULL
    bi <- NULL
    if (is.element('sm', atomic)) {
        sm <- .hedgize3(x, context=context, name=name, type='sm', allowed=hedges)
    }
    if (is.element('me', atomic)) {
        me <- .hedgize3(x, context=context, name=name, type='me', allowed=hedges)
    }
    if (is.element('bi', atomic)) {
        bi <- .hedgize3(x, context=context, name=name, type='bi', allowed=hedges)
    }
    result <- cbind.fsets(sm, me, bi)
    return(result)
}


.hedgize5 <- function(x, 
                     context,
                     name,
                     type=c('sm', 'lm', 'me', 'um', 'bi'),
                     allowed=c("ex", "ve", "ml", "ro", "ty")) {
    type <- match.arg(type)
    hedges <- .fset5[.fset5[, type], ]
    allowed <- match.arg(allowed, several.ok=TRUE)
    allowed <- c("", allowed)
    allowed <- intersect(allowed, hedges$short)

    low <- context[1]
    center <- context[2]
    big <- context[3]

    if (length(x) <= 0) {
        horizon <- numeric(0)
    } else if (type == 'sm') {
        horizon <- ifelse(low <= x & x < center,
                          (center - x) / (center - low),
                          0)
    } else if (type == 'lm') {
        lcenter <- low + (center - low) / 2
        horizon <- ifelse(low < x & x <= lcenter,
                          (x - low) / (lcenter - low),
                          ifelse(lcenter <= x & x < center,
                                 (x - center) / (lcenter - center),
                                 0))
    } else if (type == 'me') {
        horizon <- ifelse(low < x & x <= center,
                          (x - low) / (center - low),
                          ifelse(center <= x & x < big,
                                 (x - big) / (center - big),
                                 0))
    } else if (type == 'um') {
        bcenter <- center + (big - center) / 2
        horizon <- ifelse(center < x & x <= bcenter,
                          (x - center) / (bcenter - center),
                          ifelse(bcenter <= x & x < big,
                                 (x - big) / (bcenter - big),
                                 0))
    } else if (type == 'bi') {
        horizon <- ifelse(center < x & x <= big,
                          (x - center) / (big - center),
                          0)
    } else {
        horizon <- NaN
    }

    res <- NULL
    if (length(horizon) <= 0) {
        res <- matrix(0, nrow=0, ncol=length(allowed))
    } else {
        res <- sapply(allowed, function(sh) {
            # ******************************************************************************
            # ********** Pozor! tady rucne menim hedge z "me" na "ve.me"!!!!!!!!!!!!!!!!!!!!
            # ******************************************************************************
            if (type=='me' && sh=='') {
                sh <- 've'
            }
            #params <- as.vector(as.matrix((hedges[hedges$short==sh, c('p1', 'p2', 'p3')])))
            params <- as.vector(as.matrix((.fset5[.fset5$short==sh, c('p1', 'p2', 'p3')])))
            .hedge(horizon, params)
        })
    }

    if (!is.matrix(res)) {
        res <- t(as.matrix(res))
    }

    colnames(res) <- paste(.capitalize(allowed), .capitalize(type), '.', name, sep='')
    rownames(res) <- NULL

    vars <- rep(name, ncol(res))
    names(vars) <- colnames(res)

    specs <- matrix(0, nrow=nrow(.fset5), ncol=nrow(.fset5))
    specs[row(specs) > col(specs)] <- 1
    colnames(specs) <- paste(.capitalize(.fset5$short), .capitalize(type), '.', name, sep='')
    rownames(specs) <- colnames(specs)
    specs <- specs[colnames(res), colnames(res), drop=FALSE]

    return(fsets(res, vars=vars, specs=specs))
}


lcut5.numeric <- function(x, 
                          context=NULL,
                          defaultCenter=0.5,
                          atomic=c('sm', 'lm', 'me', 'um', 'bi'),
                          hedges=c("ex", "ve", "ml", "ro", "ty"),
                          name=NULL,
                          parallel=FALSE,
                          ...) {
    if (!is.vector(x) || !is.numeric(x)) {
        stop("'x' is not a numeric vector")
    }
    if (is.null(context)) {
        lo <- min(x)
        hi <- max(x)
        context <-c(lo, (hi-lo) * defaultCenter + lo, hi)
    }
    if (length(context) != 3 || context[1] != min(context) || context[3] != max(context)) {
        stop("'context' must be vector with 3 values (lo, med, hi) where lo <= med <= hi")
    }
    #if (context[1] >= context[3]) {
        #stop("'context[1]' must be lower than 'context[3]'")
    #}
    if (is.null(name)) {
        stop("If 'x' is numeric vector then 'name' must not be NULL")
    }
    if (is.null(defaultCenter) || defaultCenter < 0 || defaultCenter > 1) {
        stop("'defaultCenter' must be a number in the interval [0, 1]")
    }
    if (!is.logical(parallel) || length(parallel) != 1) {
        stop("'parallel' must be either TRUE or FALSE")
    }

    x[x < context[1]] <- context[1]
    x[x > context[3]] <- context[3]

    hedges <- match.arg(hedges, several.ok=TRUE)
    atomic <- match.arg(atomic, several.ok=TRUE)

    if (length(atomic) <= 0) {
        stop("'atomic' must not be empty")
    }

    sm <- NULL
    me <- NULL
    bi <- NULL
    if (is.element('sm', atomic)) {
        sm <- .hedgize5(x, context=context, name=name, type='sm', allowed=hedges)
    }
    if (is.element('lm', atomic)) {
        lm <- .hedgize5(x, context=context, name=name, type='lm', allowed=hedges)
    }
    if (is.element('me', atomic)) {
        me <- .hedgize5(x, context=context, name=name, type='me', allowed=hedges)
    }
    if (is.element('um', atomic)) {
        um <- .hedgize5(x, context=context, name=name, type='um', allowed=hedges)
    }
    if (is.element('bi', atomic)) {
        bi <- .hedgize5(x, context=context, name=name, type='bi', allowed=hedges)
    }
    result <- cbind.fsets(sm, lm, me, um, bi)
    return(result)
}
