###
# Generic implementation for multiple component LSS-models

### (glm/gam/m/black)boostLSS functions

mboostLSS <- function(formula, data = list(), families = GaussianLSS(),
                      control = boost_control(), weights = NULL, ...){
    cl <- match.call()
    if(is.null(cl$families))
        cl$families <- families
    fit <- mboostLSS_fit(formula = formula, data = data, families = families,
                         control = control, weights = weights, ...,
                         fun = mboost, funchar = "mboost", call = cl)
    return(fit)
}

glmboostLSS <- function(formula, data = list(), families = GaussianLSS(),
                        control = boost_control(), weights = NULL, ...){
    cl <- match.call()
    if(is.null(cl$families))
        cl$families <- families
    fit <- mboostLSS_fit(formula = formula, data = data, families = families,
                         control = control, weights = weights, ...,
                         fun = glmboost, funchar = "glmboost", call = cl)
    return(fit)
}

gamboostLSS <- function(formula, data = list(), families = GaussianLSS(),
                        control = boost_control(), weights = NULL, ...){
    cl <- match.call()
    if(is.null(cl$families))
        cl$families <- families
    fit <- mboostLSS_fit(formula = formula, data = data, families = families,
                         control = control, weights = weights, ...,
                         fun = gamboost, funchar = "gamboost", call = cl)
    return(fit)
}

blackboostLSS <- function(formula, data = list(), families = GaussianLSS(),
                          control = boost_control(), weights = NULL, ...){
    cl <- match.call()
    if(is.null(cl$families))
        cl$families <- families
    fit <- mboostLSS_fit(formula = formula, data = data, families = families,
                         control = control, weights = weights, ...,
                         fun = blackboost, funchar = "blackboost", call = cl)
    return(fit)
}


### work horse for fitting (glm/gam/m/black)boostLSS models
mboostLSS_fit <- function(formula, data = list(), families = GaussianLSS(),
                          control = boost_control(), weights = NULL,
                          fun = mboost, funchar = "mboost", call = NULL, ...){

    if (length(families) == 0)
        stop(sQuote("families"), " not specified")

    if ("offset" %in% names(list(...)))
        stop("Do not use argument ", sQuote("offset"),
             ". Please specify offsets via families")
    ### Use mu in "families" to specify offset in mu-family etc.

    if (is.list(formula)){
        if (!all(names(formula) %in% names(families)) ||
            length(unique(names(formula))) != length(names(families)))
            stop(sQuote("formula"), " can be either a formula or a named list",
                 " of formulas with same names as ",  sQuote("families"), ".")
        ynames <- sapply(formula, function(fm) as.character(fm[[2]]))
        if (length(unique(ynames)) > 1)
            warning("responses differ for the components")
        response <- lapply(formula, function(fm)
                           eval(as.expression(fm[[2]]), envir = data,
                                enclos = environment(fm)))
        #response <- eval(as.expression(formula[[1]][[2]]), envir = data,
        #                 enclos = environment(formula[[1]]))
    } else {
        response <- eval(as.expression(formula[[2]]), envir = data,
                         enclos = environment(formula))
        tmp <- vector("list", length = length(families))
        names(tmp) <- names(families)
        for (i in 1:length(tmp))
            tmp[[i]] <- formula
        formula <- tmp
    }

    mstop <- mstoparg <- control$mstop
    control$mstop <- 1
    mstop <- check(mstop, "mstop", names(families))

    nu <- control$nu
    nu <- check(nu, "nu", names(families))

    if (is.list(control$risk) || is.list(control$center) || is.list(control$trace))
        stop(sQuote("risk"),", ", sQuote("center"), " and ", sQuote("trace") ,
             " cannot be lists in ", sQuote("boost_control"))

    trace <- control$trace
    control$trace <- FALSE

    w <- weights
    if (is.null(weights)){
        if (!is.list(response)) {
            weights <- rep.int(1, NROW(response))
        } else {
            weights <- rep.int(1, NROW(response[[1]]))
        }
    }
    weights <- rescale_weights(weights)

    fit <- vector("list", length = length(families))
    names(fit) <- names(families)

    mods <- 1:length(fit)

    offset <- vector("list", length = length(mods))
    names(offset) <- names(families)
    for (j in mods){
        if (!is.list(response)) {
            response <- check_y_family(response, families[[j]])
            offset[[j]] <- families[[j]]@offset(y = response, w = weights)
        } else {
            response[[j]] <- check_y_family(response[[j]], families[[j]])
            offset[[j]] <- families[[j]]@offset(y = response[[j]], w = weights)
        }
        for (k in mods){
            for (l in mods){
                if (!is.null(offset[[l]]))
                    assign(names(offset)[l], families[[l]]@response(offset[[l]]),
                           environment(families[[k]]@ngradient))
            }
        }
    }
    for (j in mods){
        ## update value of nuisance parameters in families
         ## use response(fitted()) as this is much quicker than
         ## fitted(, type = response) as the latter uses predict()
        for (k in mods[-j]){
            if (!is.null(fit[[k]]))
                assign(names(fit)[k], families[[k]]@response(fitted(fit[[k]])),
                       environment(families[[j]]@ngradient))
        }
        ## use appropriate nu for the model
        control$nu <- nu[[j]]
        ## <FIXME> Do we need to recompute ngradient?
        fit[[j]] <- do.call(fun, list(formula[[names(families)[[j]]]],
                                      data = data, family = families[[j]],
                                      control=control, weights = w,
                                      ...))
    }
    if (trace)
        do_trace(current = 1, mstart = 0,
                 mstop = max(mstop),
                 risk = fit[[length(fit)]]$risk())

    ### set up a function for iterating boosting steps
    iBoost <- function(niter) {
        start <- sapply(fit, mstop)
        mvals <- vector("list", length(niter))
        for (j in 1:length(niter)){
            mvals[[j]] <- rep(start[j] + niter[j], max(niter))
            if (niter[j] > 0)
                mvals[[j]][1:niter[j]] <- (start[j] + 1):(start[j] + niter[j])
        }

        ENV <- lapply(mods, function(j) environment(fit[[j]]$subset))
        for (i in 1:max(niter)){
            for (j in mods){
                ## update value of nuisance parameters
                ## use response(fitted()) as this is much quicker than
                ## fitted(, type = response) as the latter uses predict()
                for (k in mods[-j])
                    assign(names(fit)[k], families[[k]]@response(fitted(fit[[k]])),
                           environment(get("ngradient", environment(fit[[j]]$subset))))
                ## update value of u, i.e. compute ngradient with new nuisance parameters

                ENV[[j]][["u"]] <- ENV[[j]][["ngradient"]](ENV[[j]][["y"]], ENV[[j]][["fit"]],
                                                           ENV[[j]][["weights"]])
                # same as:
                # evalq(u <- ngradient(y, fit, weights), environment(fit[[j]]$subset))

                ## update j-th component to "m-th" boosting step
                fit[[j]][mvals[[j]][i]]
            }
            if (trace){
                ## which is the current risk? rev() needed to get the last
                ## list element with maximum length
                whichRisk <- names(which.max(rev(lapply(lapply(fit, function(x) x$risk()), length))))
                do_trace(current = max(sapply(mvals, function(x) x[i])),
                         mstart = ifelse(firstRun, 0, max(start)),
                         mstop = ifelse(firstRun, max(niter) + 1, max(niter)),
                         risk = fit[[whichRisk]]$risk())
            }
        }
        return(TRUE)
    }

    if (any(mstop > 1)){
        ## actually go for initial mstop iterations!
        firstRun <- TRUE
        tmp <- iBoost(mstop - 1)
    }

    firstRun <- FALSE

    class(fit) <- c(paste(funchar, "LSS", sep=""), "mboostLSS")

    ### update to a new number of boosting iterations mstop
    ### i <= mstop means less iterations than current
    ### i >  mstop needs additional computations
    ### updates take place in THIS ENVIRONMENT,
    ### some models are CHANGED!
    attr(fit, "subset") <- function(i) {

        i <- check(i, "mstop", names(families))

        msf <- mstop(fit)
        niter <- i - msf
        if (all(niter == 0)) {
            ## make nothing happen with the model
            minStart <- max(msf)
        } else {
            minStart <- min(msf[niter != 0], i[niter != 0])
        }

        ## check if minStart bigger than mstop of parameters that are not touched
        #if (length(msf[niter == 0]) > 0 && minStart < min(msf[niter == 0]))
           #minStart <- min(msf[niter == 0])

        ## reduce models first (when necessary)
        if (any(msf > minStart)){
            cf <- class(fit)
            class(fit) <- "list" ## needed to use [] operator for lists

            #cat("processed parameters: ", paste(names(fit[msf > minStart]),
            #                                    collapse = ", "), "\n")

            lapply(fit[msf > minStart],
                   function(obj) obj$subset(minStart))

            ## remove additional boosting iterations from environments
            lapply(fit[msf > minStart], function(obj){
                   evalq({xselect <- xselect[1:mstop];
                          mrisk <- mrisk[1:mstop];
                          ens <- ens[1:mstop];
                          nuisance <- nuisance[1:mstop]},
                         environment(obj$subset))
               })

            class(fit) <- cf
            if (trace)
                cat("Model first reduced to mstop = ", minStart, ".\n",
                    "Now continue ...\n", sep ="")
        }

        ## now increase models (when necessary)
        if (any(i > minStart)){
            ## set negative values to 0
            ## (only applicable if some parameters do not need to be touched
            inc <- ifelse(i - minStart > 0, i - minStart, 0)
            tmp <- iBoost(inc)
        }

        mstop <<- i
    }

    ## make call in submodels nicer:
    cl <- call
    cl[[1]] <- as.name(gsub("LSS", "", cl[[1]]))
    names(cl)[names(cl) == "families"] <- "family"
    for (i in 1:length(fit)) {
        fit[[i]]$call <- cl
        ## <FIXME> This is not really nice
        fit[[i]]$call$family <- families[[i]]
    }

    attr(fit, "(weights)") <- weights  ## attach weights used for fitting

    ## update to new weights; just a fresh start
    attr(fit, "update") <- function(weights = NULL, oobweights = NULL,
                                    risk = NULL, trace = NULL, mstop = NULL) {
        if (is.null(mstop)) {
            control$mstop <- mstoparg
        } else {
            control$mstop <- mstop
        }
        if (!is.null(risk))
            control$risk <- risk
        if (!is.null(trace))
            control$trace <- trace
        ## re-use user specified offset only
        ## (since it depends on weights otherwise)
        ## this is achieved via a re-evaluation of the families argument
        mboostLSS_fit(formula = formula, data = data,
                      families = eval(call[["families"]]), weights = weights,
                      control = control, fun = fun, funchar = funchar,
                      call = call, oobweights = oobweights)
    }
    attr(fit, "control") <- control
    attr(fit, "call") <- call
    attr(fit, "data") <- data
    attr(fit, "families") <- families
    return(fit)
}
