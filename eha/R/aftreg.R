aftreg <- function (formula = formula(data),
                    data = parent.frame(),
                    na.action = getOption("na.action"),
                    dist = "weibull",
                    init,
                    shape = 0,
                    id,
                    ##param = c("default", "survreg", "canonical"),
                    param = c("lifeAcc", "lifeExp"),
                    control = list(eps = 1e-8, maxiter = 20, trace = FALSE),
                    singular.ok = TRUE,
                    model = FALSE,
                    x = FALSE,
                    y = TRUE)
{
    param <- param[1]
    if (param == "survreg"){
        param <- "lifeExp" # backwards compability
        warning("'survreg' is a deprecated argument value")
    }else if (param == "canonical"){
        param <- "lifeExp"
        warning("'canonical' is a deprecated argument value")
    }else{
        if (!(param %in% c("lifeAcc", "lifeExp"))){
            stop(paste(param, "is not a valid parametrization."))
        }
    }
    ## if (dist == "gompertz") shape <- 1
    pfixed <- any(shape > 0)
    call <- match.call()
    m <- match.call(expand.dots = FALSE)

    temp <- c("", "formula", "data", "id", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)] # m is a call

    special <- "strata"
    Terms <- if (missing(data))
        terms(formula, special)
    else terms(formula, special, data = data) # Terms is a 'terms' 'formula'

    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    ##return(m)
        Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv"))
        stop("Response must be a survival object")
    if (missing(id)) id <- 1:nrow(Y)
    else id <- model.extract(m, "id")
    ##else id <- m$"(id)" # This does not work; leave it for the time being...
    
    ##weights <- model.extract(m, "weights") # No weights (as yet...)
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
    ##return(X)
    assign <- lapply(attrassign(X, newTerms)[-1], function(x) x - 1)
    
    X <- X[, -1, drop = FALSE]  ##OBS!!!! No Intercept!
    
    ncov <- NCOL(X)
    
#########################################
    
    if (ncov){
        if (length(dropx)){
            covars <- names(m)[-c(1, (dropx + 1))]
        }else{
            covars <- names(m)[-1]
        }
        
        isF <- logical(length(covars))
        for (i in 1:length(covars)){
            if (length(dropx)){
                if (is.logical(m[, -(dropx + 1)][, (i + 1)])){
                    m[, -(dropx + 1)][, (i + 1)] <-
                        as.factor(m[, -(dropx + 1)][, (i + 1)])
                }
                isF[i] <- is.factor(m[, -(dropx + 1)][, (i + 1)])## ||
                ##is.logical(m[, -(dropx + 1)][, (i + 1)]) )
            }else{
                if (is.logical(m[, (i + 1)])){
                    m[, (i + 1)] <- as.factor(m[, (i + 1)])
                }
                isF[i] <- is.factor(m[, (i + 1)]) ##||
                ##is.logical(m[, (i + 1)]) )
            }
        }
        
        if (ant.fak <- sum(isF)){
            levels <- list()
            index <- 0
            for ( i in 1:length(covars) ){
                if (isF[i]){
                    index <- index + 1
                    if (length(dropx)){
                        levels[[i]] <- levels(m[, -(dropx + 1)][, (i + 1)])
                    }else{
                        levels[[i]] <- levels(m[, (i + 1)])
                    }
                }else{
                    levels[[i]] <- NULL
                }
            }
        }else{
            levels <- NULL
        }
    }
    
##########################################
    type <- attr(Y, "type")
    if (type != "right" && type != "counting")
        stop(paste("This model doesn't support \"", type, "\" survival data",
                   sep = ""))
    
    if (NCOL(Y) == 2){
        Y <- cbind(numeric(NROW(Y)), Y)
    }
    
    n.events <- sum(Y[, 3] != 0)
    if (n.events == 0) stop("No events; no sense in continuing!")
    if (missing(init))
        init <- NULL
    
    if (is.list(control)){
        if (is.null(control$eps)) control$eps <- 1e-8
        if (is.null(control$maxiter)) control$maxiter <- 10
        if (is.null(control$trace)) control$trace <- FALSE
    }else{
        stop("control must be a list")
    }
    
    ##cat("\nEntering aftreg.fit ...")
    fit <- aftreg.fit(X,
                      Y,
                      dist,
                      param,
                      strats,
                      offset,
                      init,
                      shape,
                      id,
                      control,
                      pfixed)
    ##cat("and back!\n\n")
    if (!is.null(fit$overlap)) return(fit$overlap)
    
    if (ncov){
        fit$linear.predictors <- offset + X %*% fit$coefficients[1:ncov]
        fit$means <- apply(X, 2, mean)
    }else{
        fit$linear.predictors <- offset
        fit$means <- NULL
    }
    ##score <- exp(lp)



    if (!fit$fail){
        fit$fail <- NULL
    }else{
        out <- paste("Singular hessian; suspicious variable No. ",
                     as.character(fit$fail), ":\n",
                     names(coefficients)[fit$fail], " = ",
                     as.character(fit$value),
                     "\nTry running with fixed shape", sep = "")
        stop(out)
    }


    fit$convergence <- as.logical(fit$conver)
    fit$conver <- NULL ## Ugly!

###########################################################################
    ## Crap dealt with ......

    if (is.character(fit)) {
        fit <- list(fail = fit)
        class(fit) <- "mlreg"
    }
    else if (is.null(fit$fail)){
        if (!is.null(fit$coef) && any(is.na(fit$coef))) {
            vars <- (1:length(fit$coef))[is.na(fit$coef)]
            msg <- paste("X matrix deemed to be singular; variable",
                         paste(vars, collapse = " "))
            if (singular.ok)
              warning(msg)
            else stop(msg)
        }
        fit$n <- nrow(Y)
        fit$terms <- Terms
        fit$assign <- assign
        if (FALSE){ ## Out-commented...(why?)
            if (length(fit$coef) && is.null(fit$wald.test)) {
                nabeta <- !is.na(fit$coef)
                if (is.null(init))
                  temp <- fit$coef[nabeta]
                else temp <- (fit$coef - init)[nabeta]
                ##fit$wald.test <-
                  ##survival:::coxph.wtest(fit$var[nabeta, nabeta],
                    ##                     temp, control$toler.chol)$test
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
        if (y)
            fit$y <- Y
    }
    ##if (!is.null(weights) && any(weights != 1))
    ##    fit$weights <- weights

##########################################
    s.wght <- (Y[, 2] - Y[, 1])## * weights
    fit$ttr <- sum(s.wght)
    if (ncov){
        fit$isF <- isF
        fit$covars <- covars
        fit$w.means <- list()
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
        fit$means <- colMeans(X)

    }
    ##cat("\nHere we go !!!!!!!!!!!!!!\n\n")
##########################################
    fit$ttr <- sum(s.wght)
    ##names(fit$coefficients) <- coef.names
    fit$levels <- levels
    fit$formula <- formula(Terms)
    fit$call <- call
    fit$dist <- dist
    fit$events <- n.events

    class(fit) <- c("aftreg", "phreg")
    fit$param <- param # New in 2.1-1:
    if (dist == "gompertz"){
        baselineMean <- numeric(fit$n.strata)
        for (j in 1:fit$n.strata){
            scale <- exp(fit$coef[ncov + 2 * j - 1])
            shape <- exp(fit$coef[ncov + 2 * j])
            ## Simulation!
            baselineMean[j] <- mean(rgompertz(100000, param = "canonical",
                                              scale = scale, shape = shape))
        }
        fit$baselineMean <- baselineMean
    }else{ # To be filled for other dists!
        fit$baselineMean <- NULL
    }

    ##
    fit$pfixed <- pfixed
    fit
}
