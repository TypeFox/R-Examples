coxreg <- function (formula = formula(data),
                    data = parent.frame(),
                    weights,
                    subset,
                    t.offset,
                    na.action = getOption("na.action"),
                    init = NULL,
                    method = c("efron", "breslow", "mppl", "ml"),
                    control = list(eps = 1e-8, maxiter = 25, trace = FALSE),
                    singular.ok = TRUE,
                    model = FALSE,
                    center = TRUE,
                    x = FALSE,
                    y = TRUE,
                    hazards = TRUE,
                    boot = FALSE,
                    efrac = 0,
                    geometric = FALSE,
                    rs = NULL,
                    frailty = NULL,
                    max.survs = NULL)
{

    meth <- method[1]
    cox.ph <- (missing(t.offset) &&
               (meth %in% c("breslow", "efron")) &&
               is.null(rs) &&
               is.null(max.survs) &&
               (!boot) &&
               (efrac == 0) &&
               is.null(frailty) &&
               (!geometric)
               )
    
    if (!is.null(frailty))
        stop("Frailty not implemented (yet). Try the 'coxme' package")
    method <- match.arg(method)
    
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data", "weights", "subset", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]

    special <- "strata"
    Terms <- if (missing(data))
        terms(formula, special)
    else
        terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())

    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv"))
        stop("Response must be a survival object")

    if (is.null(max.survs)) max.survs <- NROW(Y)
    if (missing(weights)) weights <- rep(1, NROW(Y))
    else weights <- model.extract(m, "weights")
    cox.ph <- cox.ph && (length(weights) == NROW(Y))
    
    if (missing(t.offset)) t.offset <- NULL
    ##
    offset <- attr(Terms, "offset")
    tt <- length(offset)
    offset <- if (tt == 0)
        rep(0, nrow(Y))
    else if (tt == 1)
        m[[offset]]
    else {
        ff <- m[[offset[1]]]
        for (i in 2:tt) ff <- ff + m[[offset[i]]]
        ff
    }
    attr(Terms, "intercept") <- 1
    strats <- attr(Terms, "specials")$strata
    dropx <- NULL
        
    if (length(strats)) {
        temp <- untangle.specials(Terms, "strata", 1)
        dropx <- c(dropx, temp$terms)
        if (length(temp$vars) == 1)
            strata.keep <- m[[temp$vars]]
        else strata.keep <- strata(m[, temp$vars], shortlabel = TRUE)
        strats <- as.numeric(strata.keep)
    }
    if (length(dropx))
        newTerms <- Terms[-dropx]
    else newTerms <- Terms
    X <- model.matrix(newTerms, m)
    assign <- lapply(attrassign(X, newTerms)[-1], function(x) x -
                     1)
    X <- X[, -1, drop = FALSE]
    ##}
#########################################
        
    if (length(dropx)){
        covars <- names(m)[-c(1, (dropx + 1))]
    }else{
        covars <- names(m)[-1]
    }

    isI <- logical(NCOL(X)) # Added Jan 2014; 2.4-0.
    isF <- logical(length(covars))
    isO <- logical(length(covars))
    if (length(covars)){
        for (i in 1:length(covars)){
            if (length(dropx)){
                if (is.logical(m[, -(dropx + 1)][, (i + 1)])){
                    m[, -(dropx + 1)][, (i + 1)] <-
                        as.factor(m[, -(dropx + 1)][, (i + 1)])
                }
                isF[i] <- is.factor(m[, -(dropx + 1)][, (i + 1)])## ||
                isO[i] <- is.ordered(m[, -(dropx + 1)][, (i + 1)])## ||
                ##is.logical(m[, -(dropx + 1)][, (i + 1)]) )
            }else{
                if (is.logical(m[, (i + 1)])){
                    m[, (i + 1)] <- as.factor(m[, (i + 1)])
                }
                isF[i] <- is.factor(m[, (i + 1)]) ##||
                isO[i] <- is.ordered(m[, (i + 1)]) ##||
                ## is.logical(m[, (i + 1)]) )
            }
        }
    }
    
    if (any(isF)){
        levels <- list()
        index <- 0
        for ( i in 1:length(covars) ){
            if (isF[i]){
                index <- index + 1
                if (length(dropx)){
                    ll <- levels(m[, -(dropx + 1)][, (i + 1)])
                    if (isO[i]){
                        levels[[i]] <- paste("order", seq_along(ll) - 1)
                    }else{
                        levels[[i]] <- ll
                    }
                }else{
                    ll <- levels(m[, (i + 1)])
                    if (isO[i]){
                        levels[[i]] <- paste("order", seq_along(ll) - 1)
                    }else{
                        levels[[i]] <- ll
                    }
                }
            }else{
                levels[[i]] <- NULL
            }
        }
    }else{
        levels <- NULL
    }

   if (any(isF)){ ## New; get isI: (Jan 2014; 2.4-0)
        indx <- 0
        for (i in 1:length(covars)){
            indx <- indx + 1
            if (isF[i]){
                isI[indx] <- TRUE
                if (length(levels[[i]]) >= 3){
                    for (j in 3:length(levels[[i]])){
                        indx <- indx + 1
                        isI[indx] <- TRUE
                    }
                }
            }
        }
    }

    n.events <- sum(Y[, NCOL(Y)] != 0)
    if (n.events == 0) stop("No events; no sense in continuing!")
        
    ##########################################

    ## Fixed now? if (FALSE){      ## This has to be fixed in the future!!
    if (NCOL(X) == 0){ # No covariates; special treatment!
        if (is.null(strats)) stratum <- rep(1, NROW(Y))
        else stratum <- strats
        type <- attr(Y, "type")
        control$iter.max <- 0
        control$toler.chol <- .Machine$double.eps^0.75
        control$toler.inf <- sqrt(control$eps)
        control$outer.max <- 10
        X <- matrix(0, nrow = NROW(Y), ncol = 1)
        init <- 0
        if (type == "counting"){
            fit <- survival::agreg.fit(X, Y, stratum, offset, init,
                                        control, weights = weights,
                                        method = method, row.names(m))
        }else{
            fit <- survival::coxph.fit(X, Y, stratum, offset, init,
                                        control, weights = weights,
                                        method = method, row.names(m))
        }
        fit$nullModel <- TRUE
        if (hazards){
            scores <- exp(offset)
            hazards <- getHaz(Y, stratum, scores)
            class(hazards) <- "hazdata"
            fit$hazards <- hazards
        }
        fit$call <- call
        fit$n.events <- n.events
        fit$n <- NROW(Y)
        fit$y <- Y
        class(fit) <- c("coxreg", "coxph")
        fit$means <- 0
        return(fit)
    
    }else if (cox.ph){
        type <- attr(Y, "type")
        control$iter.max <- control$maxiter
        control$toler.chol <- .Machine$double.eps^0.75
        control$toler.inf <- sqrt(control$eps)
        control$outer.max <- 10
        if (is.null(strats)) stratum <- rep(1, NROW(Y))
        else stratum <- strats
        if (type == "counting"){
            fit <- survival::agreg.fit(X, Y, stratum, offset, init,
                                        control, weights = weights,
                                        method = method, row.names(m))
        }else{
            fit <- survival::coxph.fit(X, Y, stratum, offset, init,
                                        control, weights = weights,
                                        method = method, row.names(m))
        }
        fit$nullModel <- FALSE
        ## get hazards
        ## New in 2.4-0: covariates are centered; indicators not!
        ## If center == FALSE, X.means are "added back" before call to
        ## getHaz!!

        if (center){
            X.means <- colMeans(X)
            for (i in seq_len(NCOL(X))){
                if (isI[i]) X.means[i] <- 0
            }
            scores <- exp(offset + X %*% fit$coefficients -
                          sum(X.means * fit$coefficients))
        }else{
            X.means <- numeric(NCOL(X))
            scores <- exp(offset + X %*% fit$coefficients)
        }

        if (hazards){
            hazards <- getHaz(Y, stratum, scores)
            class(hazards) <- "hazdata"
            fit$hazards <- hazards
        }else{
            fit$hazards <- NULL
        }
            ##rs <- risksets(Y, strats)
        ##hazard <- .Fortran("gethaz",
        ##                   as.integer(NROW(Y)),  # 'nn'
        ##                   as.integer(length(rs$antrs)), # 'ns' 
        ##                   as.integer(rs$antrs), # 'antrs'
        ##                   as.integer(rs$size), # 'size'
        ##                   as.integer(rs$n.events), # 'nevents'
        ##                   as.integer(length(rs$riskset)), # 'totsize'
        ##                   as.integer(rs$riskset), # 'riskset'
        ##                   as.double(exp(fit$linear.predictors)), # 'score'
        ##                   as.integer(sum(rs$antrs)), # 'totrs'
        ##                   hazard = double(sum(rs$antrs)), # 'hazard' (return)
        ##                   DUP = FALSE,
        ##                   PACKAGE = "eha")$hazard
        ## Put it on:
        ##haz.mean <- fit$hazard::: At means of covariates:
        ##if (!is.null(fit$coefficients))
          ##  hazard <- 1 - (1 - hazard)^exp(fit$means * fit$coefficients)
        ##hazards <- list()
        ##stopp <- cumsum(rs$antrs)
        ##startt <- c(1, 1 + stopp[-length(rs$antrs)])
        ##for (i in 1:length(rs$antrs)){
          ##  hazards[[i]] <- cbind(rs$risktimes[startt[i]:stopp[i]],
            ##                      hazard[startt[i]:stopp[i]])
        ##}
        ##fit$hazards <- hazard
    }else{ # if (!cox.ph)
        if (NCOL(Y) == 2){
            Y <- cbind(numeric(NROW(Y)), Y)
            attr(Y, "type") <- "counting"
        }
        
        ##return(Y)
        type <- attr(Y, "type")
        if (type != "right" && type != "counting")
            stop(paste("Cox model doesn't support \"", type, "\" survival data",
                       sep = ""))
        
        if ((!is.null(init)) && (length(init) != NCOL(X)))
            stop("Wrong length of 'init'")
        
        
        if (is.list(control)){
            if (is.null(control$eps)) control$eps <- 1e-8
            if (is.null(control$maxiter)) control$maxiter <- 10
            if (is.null(control$trace)) control$trace <- FALSE
        }else{
            stop("control must be a list")
        }
        
### New start for cox.ph (not any more!) ##################################
        if (geometric){
            method <- "ml"
            fit <- geome.fit(X,
                             Y,
                             rs,
                             strats,
                             offset,
                             init,
                             max.survs,
                             method,
                             boot,
                             control)
            fit$nullModel <- FALSE
        }else{
            fit <- coxreg.fit(X,
                              Y,
                              rs,
                              weights,
                              t.offset,
                              strats,
                              offset,
                              init,
                              max.survs,
                              method,
                              center,
                              boot,
                              efrac,
                              calc.hazards = FALSE, ## Changed 2.4-0
                              calc.martres = TRUE,
                              control,
                              verbose = TRUE)
            ## get hazards
            fit$nullModel <- FALSE
            if (center){
                X.means <- colMeans(X)
                for (i in seq_len(NCOL(X))){
                    if (isI[i]) X.means[i] <- 0
                }
                scores <- exp(offset + X %*% fit$coefficients -
                              sum(X.means * fit$coefficients))
            }else{
                X.means <- numeric(NCOL(X))
                scores <- exp(offset + X %*% fit$coefficients)
            }

            if (is.null(strats)){
                stratum <- rep(1, NROW(Y))
            }else{
                stratum <- strats
                fit$stratum <- strats
            }
            if (hazards){
                hazards <- getHaz(Y, stratum, scores)
                class(hazards) <- "hazdata"
                fit$hazards <- hazards
            }else{
                fit$hazards <- NULL
            }
        }

        fit$convergence <- as.logical(fit$conver)
        fit$conver <- NULL ## Ugly!
        fit$f.convergence <- as.logical(fit$f.conver)
        fit$f.conver <- NULL
    }
###########################################################################
## Crap dealt with ......

    if (is.character(fit)) {
        fit <- list(fail = fit)
        class(fit) <- "coxreg"
    }else if ((!cox.ph) && !fit$fail){
        if (length(fit$coef) && any(is.na(fit$coef))) {
            vars <- (1:length(fit$coef))[is.na(fit$coef)]
            msg <- paste("X matrix deemed to be singular; variable",
                         paste(vars, collapse = " "))
            if (singular.ok)
                warning(msg)
            else stop(msg)
        }
        fit$n <- nrow(Y)
        class(fit) <- fit$method
        fit$terms <- Terms
        fit$assign <- assign
        if (FALSE){ ## Out-commented
            if (length(fit$coef) && is.null(fit$wald.test)) {
                nabeta <- !is.na(fit$coef)
                if (is.null(init))
                    temp <- fit$coef[nabeta]
                else temp <- (fit$coef - init)[nabeta]
                ##fit$wald.test <-
                  ##  survival:::coxph.wtest(fit$var[nabeta, nabeta],
                    ##                       temp, control$toler.chol)$test
            }
        }
        na.action <- attr(m, "na.action")
        if (length(na.action))
            fit$na.action <- na.action
        if (model)
            fit$model <- m
        if (x) {
            fit$x <- X
            if (length(strats))
                fit$strata <- strata.keep
        }
    }
    fit$stratum <- strats
    if (y)
        fit$y <- Y
    if (x)
        fit$x <- X
    ##if (!is.null(weights) && any(weights != 1))
    ##    fit$weights <- weights

    ##########################################

    fit$isI <- isI
    fit$isF <- isF
    fit$isO <- isO
    fit$covars <- covars
    if (NCOL(Y) == 3){
        s.wght <- (Y[, 2] - Y[, 1])## * weights
    }else{
        s.wght <- Y[, 1]
    }
    fit$ttr <- sum(s.wght)
    fit$w.means <- list()
    if (length(fit$covars)){
        for (i in 1:length(fit$covars)){
            nam <- fit$covars[i]
            col.m <- which(nam == names(m))
            if (isF[i]){
                n.lev <- length(levels[[i]])
                fit$w.means[[i]] <- numeric(n.lev)
                for (j in 1:n.lev){
                    who <- m[, col.m] == levels[[i]][j]
                    fit$w.means[[i]][j] <-
                      sum( s.wght[who] ) / fit$ttr ## * 100, if in per cent
                }
            }else{
                fit$w.means[[i]] <- sum(s.wght * m[, col.m]) / fit$ttr
            }
        }
    }
    ##########################################
    fit$nullModel <- FALSE
    fit$levels <- levels
    fit$formula <- formula(Terms)
    fit$terms <- Terms
    fit$call <- call
    fit$events <- n.events
    fit$center <- center
    if (length(fit$coefficients)){
        names(fit$coefficients) <- colnames(X)
        fit$means <- X.means
    }else{
        fit$means <- numeric(0)
    }
    if (length(strats)){
        fit$strata <- levels(as.factor(strata.keep)) ## New 2.2-6
    }
    fit$method <- method
    fit$n <- NROW(Y)
    fit$df <- length(fit$coefficients)
    class(fit) <- c("coxreg", "coxph") # Not Removed "coxph"; cox.zph!
    ##class(fit) <- "coxreg"
    fit
}
