`cusp` <-
function (formula, alpha, beta, data, weights, offset, ..., control = glm.control(), 
    method='cusp.fit', optim.method='L-BFGS-B', model = TRUE, contrasts = NULL) 
{
    call <- match.call()
    if (missing(data)) 
        data <- environment(formula)
    # model matrix y
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    if(length(mf[[2]])>2) { mf[[2]][[2]] <- NULL }
    mf <- eval(mf, envir = parent.frame())
    Y <- model.response(mf, "any")
    mt <- attr(mf,'terms')
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0)
    
    # model matrix alpha
    formula.alpha <- alpha
    mfa = match.call(expand.dots=FALSE)
    m <- match(c("alpha","data","subset","weights","na.action","offset"), names(mfa), 0)
    mfa <- mfa[c(1,m)]
    mfa$drop.unused.levels <- TRUE
    mfa[[1]] <- as.name("model.frame")
    mfa[[2]] <- update(alpha, paste(attr(terms(formula),'term.labels')[1],"~ ."))
    names(mfa) <- c('','formula',names(mfa)[-(1:2)])
    mfa <- eval(mfa, envir = parent.frame()) # is nodig als data argument niet is mee gegeven lijkt me...
    mta <- attr(mfa,'terms')
    X.alpha <- if (!is.empty.model(mta)) 
        model.matrix(mta, mfa, contrasts)
 
    # model matrix beta
    formula.beta <- beta
    mfb = match.call(expand.dots=FALSE)
    m <- match(c("beta","data","subset","weights","na.action","offset"), names(mfb), 0)
    mfb <- mfb[c(1,m)]
    mfb$drop.unused.levels <- TRUE
    mfb[[1]] <- as.name("model.frame")
    mfb[[2]] <- update(beta, paste(attr(terms(formula),'term.labels')[1],"~ ."))
    names(mfb) <- c('','formula',names(mfb)[-(1:2)])
    mfb <- eval(mfb, envir = parent.frame()) # is nodig als data argument niet is mee gegeven lijkt me...
    mtb <- attr(mfb,'terms')
    X.beta <- if (!is.empty.model(mtb)) 
        model.matrix(mtb, mfb, contrasts)
    
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1) 
            offset <- rep(offset, NROW(Y))
        else if (length(offset) != NROW(Y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
    }
    fit <- cusp.fit(y = Y <- X, x.alpha = X.alpha, x.beta = X.beta, ..., method=optim.method)
    if (model){
        fit$model <- list(y = mf, alpha = mfa, beta  = mfb)
    }
    fit$na.action <- attr(mf, "na.action")
    fit$x <- list(X.y = X, X.alpha = X.alpha, X.beta = X.beta); 
        
    fit <- c(fit, list(call = call, formula = list(y=formula, alpha=formula.alpha, 
        beta=formula.beta), terms = list(y = mt, alpha = mta, beta = mtb),
        OK = all(eigen(fit$hess, symmetric=TRUE, only.values=TRUE)$values>0),
        data = data, offset = offset, control = control, method = method,
        contrasts = list(y=attr(X, "contrasts"), alpha=attr(X.beta, "contrasts"),
        beta=attr(X.beta, "contrasts")), xlevels = list(y=.getXlevels(mt, mf),
        alpha=.getXlevels(mta, mfa), beta=.getXlevels(mtb, mfb))))
    class(fit) <- c("cusp","glm","lm")
    fit
}

