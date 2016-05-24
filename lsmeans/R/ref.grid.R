# Reference grid code


# Change to cov.reduce specification: can be...
#     a function: is applied to all covariates
#     named list of functions: applied to those covariates (else mean is used)
#     TRUE - same as mean
#     FALSE - same as function(x) sort(unique(x))

ref.grid <- function(object, at, cov.reduce = mean, mult.name, mult.levs, 
                     options = get.lsm.option("ref.grid"), data, type, ...) {
    # recover the data
    if (missing(data)) {
        data = try(recover.data (object, data = NULL, ...))
        if (inherits(data, "try-error"))
            stop("Perhaps a 'data' or 'params' argument is needed")
    }
    else # attach needed attributes to given data
        data = recover.data(object, data = data, ...)
    
    if(is.character(data)) # 'data' is in fact an error message
        stop(data)
        
    
    trms = attr(data, "terms")
    
    # find out if any variables are coerced to factors
    ### OLD VERSION: anm = all.names(attr(data, "terms"))    
    ###              coerced = anm[1 + grep("factor|ordered", anm)]
    coerced = .find.coerced(trms, data)
    
    # convenience function
    sort.unique = function(x) sort(unique(x))
    
    # Ensure cov.reduce is a function or list thereof
    dep.x = list() # list of formulas to fit later
    fix.cr = function(cvr) {
        # cvr is TRUE or FALSE
        if(is.logical(cvr)) 
            if(cvr[1]) cvr = mean
        else              cvr = sort.unique
        else if (inherits(cvr, "formula")) {
            if (length(cvr) < 3)
                stop("Formulas in 'cov.reduce' must be two-sided")
            lhs = .all.vars(cvr)[1]
            dep.x[[lhs]] <<- cvr
            cvr = mean 
        }
        else if (!inherits(cvr, c("function","list")))
            stop("Invalid 'cov.reduce' argument")
        cvr
    }
    
    # IMPORTANT: following stmts may also affect x.dep
    if (is.list(cov.reduce))
        cov.reduce = lapply(cov.reduce, fix.cr)
    else
        cov.reduce = fix.cr(cov.reduce)
    
    # zap any formulas that are also in 'at'
    if (!missing(at))
        for (xnm in names(at)) dep.x[[xnm]] = NULL
    
    
    # local cov.reduce function that works with function or named list
    cr = function(x, nm) {
        if (is.function(cov.reduce))
            cov.reduce(x)
        else if (!is.null(cov.reduce[[nm]]))
            cov.reduce[[nm]](x)
        else
            mean(x)
    }
    
    # initialize empty lists
    ref.levels = matlevs = xlev = list()
    
    for (nm in attr(data, "responses")) {
        y = data[[nm]]
        if (is.matrix(y))
            matlevs[[nm]] = apply(y, 2, mean)
        else
            ref.levels[[nm]] = mean(y)
    }
    
    for (nm in attr(data, "predictors")) {
        x = data[[nm]]
        
        # Save the original levels of factors, no matter what
        if (is.factor(x))
            xlev[[nm]] = levels(factor(x))
            # (applying factor drops any unused levels)
    
        # Now go thru and find reference levels...
        # mentioned in 'at' list but not coerced
        if (!(nm %in% coerced) && !missing(at) && !is.null(at[[nm]]))
            ref.levels[[nm]] = at[[nm]]
        # factors not in 'at'
        else if (is.factor(x))
            ref.levels[[nm]] = levels(factor(x))
        else if (is.character(x))
            ref.levels[[nm]] = sort.unique(x)
        # matrices
        else if (is.matrix(x)) {
            # Matrices -- reduce columns thereof, but don't add to baselevs
            matlevs[[nm]] = apply(x, 2, cr, nm)
            # if cov.reduce returns a vector, average its columns
            if (is.matrix(matlevs[[nm]]))
                matlevs[[nm]] = apply(matlevs[[nm]], 2, mean)
        }
        # covariate coerced, or not mentioned in 'at'
        else {
            # single numeric pred but coerced to a factor - use unique values
            # even if in 'at' list. We'll fix this up later
            if (nm %in% coerced)            
                ref.levels[[nm]] = sort.unique(x)
            
            # Ordinary covariates - summarize
            else 
                ref.levels[[nm]] = cr(x, nm)
        }
    }
    
    # Now create the reference grid
    grid = do.call(expand.grid, ref.levels)
    
    # add any matrices
    for (nm in names(matlevs))
        grid[[nm]] = matrix(rep(matlevs[[nm]], each=nrow(grid)), nrow=nrow(grid))

    # resolve any covariate formulas
    for (xnm in names(dep.x)) {
        if (!all(.all.vars(dep.x[[xnm]]) %in% names(grid)))
            stop("Formulas in 'cov.reduce' must predict covariates actually in the model")
        xmod = lm(dep.x[[xnm]], data = data)
        grid[[xnm]] = predict(xmod, newdata = grid)
        ref.levels[[xnm]] = NULL
    }
    
    basis = lsm.basis(object, trms, xlev, grid, ...)
    
    misc = basis$misc
    
    form = attr(data, "call")$formula
    if (is.null(misc$tran) && (length(form) > 2)) { # No link fcn, but response may be transformed
        lhs = form[-3] ####form[[2]]
        tran = setdiff(.all.vars(lhs, functions = TRUE), c(.all.vars(lhs), "~", "cbind"))
        if(length(tran) > 0) {
            tran = paste(tran, collapse = ".")  
            # length > 1: Almost certainly unsupported, but facilitates a more informative error message
            
            # Look for a multiplier, e.g. 2*sqrt(y)
            tst = strsplit(strsplit(as.character(form[2]), "\\(")[[1]][1], "\\*")[[1]]
            if(length(tst) > 1) {
                mul = suppressWarnings(as.numeric(tst[1]))
                if(!is.na(mul))
                    misc$tran.mult = mul
                tran = gsub("\\*\\.", "", tran)
            }
            if (tran == "linkfun")
                tran = as.list(environment(trms))
            misc$tran = tran
            misc$inv.lbl = "response"
        }
    }
    
    # Take care of multivariate response
    multresp = character(0) ### ??? was list()
    ylevs = misc$ylevs
    if(!is.null(ylevs)) { # have a multivariate situation
       if (missing(mult.levs)) {
            if (missing(mult.name))
                mult.name = names(ylevs)[1]
            ref.levels[[mult.name]] = ylevs[[1]]
            multresp = mult.name
            MF = data.frame(ylevs)
            names(MF) = mult.name
        }
        else {
            k = prod(sapply(mult.levs, length))
            if (k != length(ylevs[[1]])) 
                stop("supplied 'mult.levs' is of different length than that of multivariate response")
            for (nm in names(mult.levs))
                ref.levels[[nm]] = mult.levs[[nm]]
            multresp = names(mult.levs)
            MF = do.call("expand.grid", mult.levs)
        }
        ###grid = do.call("expand.grid", ref.levels)
        grid = merge(grid, MF)
        # add any matrices
        for (nm in names(matlevs))
            grid[[nm]] = matrix(rep(matlevs[[nm]], each=nrow(grid)), nrow=nrow(grid))
    }

# Here's a complication. If a numeric predictor was coerced to a factor, we had to
# include all its levels in the reference grid, even if altered in 'at'
# Moreover, whatever levels are in 'at' must be a subset of the unique values
# So we now need to subset the rows of the grid and linfct based on 'at'

    problems = if (!missing(at)) 
        intersect(c(multresp, coerced), names(at)) 
    else character(0)
    if (length(problems > 0)) {
        incl.flags = rep(TRUE, nrow(grid))
        for (nm in problems) {
            if (is.numeric(ref.levels[[nm]])) {
                at[[nm]] = round(at[[nm]], 3)
                ref.levels[[nm]] = round(ref.levels[[nm]], 3)
            }
            # get only "legal" levels
            at[[nm]] = at[[nm]][at[[nm]] %in% ref.levels[[nm]]]
            # Now which of those are left out?
            excl = setdiff(ref.levels[[nm]], at[[nm]])
            for (x in excl)
                incl.flags[grid[[nm]] == x] = FALSE
            ref.levels[[nm]] = at[[nm]]
        }
        if (!any(incl.flags))
            stop("Reference grid is empty due to mismatched levels in 'at'")
        grid = grid[incl.flags, , drop=FALSE]
        basis$X = basis$X[incl.flags, , drop=FALSE]
    }

    # Any offsets??? (misc$offset.mult might specify removing or reversing the offset)
    if(!is.null(attr(trms,"offset"))) {
        om = 1
        if (!is.null(misc$offset.mult))
            om = misc$offset.mult
        if (any(om != 0))
            grid[[".offset."]] = om * .get.offset(trms, grid)
    }

    ### --- Determine weights for each grid point --- (added ver.2.11), updated ver.2.14 to include weights
    if (is.null(data[["(weights)"]]))
        data[["(weights)"]] = 1
    nms = union(names(xlev), coerced) # only factors, no covariates or mult.resp
    # originally, I used 'plyr::count', but there are probs when there is a 'freq' variable
    id = plyr::id(data[, nms, drop = FALSE], drop = TRUE)
    uid = !duplicated(id)
    key = do.call(paste, data[uid, nms, drop = FALSE])
    key = key[order(id[uid])]
    #frq = tabulate(id, attr(id, "n"))
    tgt = do.call(paste, grid[, nms, drop = FALSE])
    wgt = rep(0, nrow(grid))
    for (i in seq_along(key))
        wgt[tgt == key[i]] = sum(data[["(weights)"]][id==i])
    grid[[".wgt."]] = wgt

    misc$ylevs = NULL # No longer needed
    misc$estName = "prediction"
    misc$estType = "prediction"
    misc$infer = c(FALSE,FALSE)
    misc$level = .95
    misc$adjust = "none"
    misc$famSize = nrow(grid)
    misc$avgd.over = character(0)

    post.beta = basis$post.beta
    if (is.null(post.beta))
        post.beta = matrix(NA)
    
    result = new("ref.grid",
         model.info = list(call = attr(data,"call"), terms = trms, xlev = xlev),
         roles = list(predictors = attr(data, "predictors"), 
                      responses = attr(data, "responses"), 
                      multresp = multresp),
         grid = grid, levels = ref.levels, matlevs = matlevs,
         linfct = basis$X, bhat = basis$bhat, nbasis = basis$nbasis, V = basis$V,
         dffun = basis$dffun, dfargs = basis$dfargs, 
         misc = misc, post.beta = post.beta)
        
    if (!missing(type)) {
        if (is.null(options)) options = list()
        options$predict.type = type
    }

    if(!is.null(options)) {
        options$object = result
        result = do.call("update.ref.grid", options)
    }

    if(!is.null(hook <- misc$postGridHook)) {
        if (is.character(hook))
            hook = get(hook)
        result@misc$postGridHook = NULL
        result = hook(result)
    }
    
    .save.ref.grid(result)
    result
}


#### End of ref.grid ------------------------------------------

# local utility to identify ref.grid that is not an lsmobj
.is.true.ref.grid = function(object) {
    is(object, "ref.grid") && !is(object, "lsmobj")
}

# local utility to save each newly constructed ref.grid, if enabled
# Goes into global environment unless .Last.ref.grid is found further up
.save.ref.grid = function(object) {
    if(get.lsm.option("save.ref.grid") && .is.true.ref.grid(object))
        assign(".Last.ref.grid", object, inherits = TRUE)
}



# This function figures out which covariates in a model 
# have been coerced to factors. Does NOT rely on the names of
# functions like 'factor' or 'interaction' as we use actual results
.find.coerced = function(trms, data) {
    isfac = sapply(data, function(x) inherits(x, "factor"))
    
    # Character vectors of factors and covariates in the data...
    facs.d = names(data)[isfac]
    covs.d = names(data)[!isfac]
    
    lbls = attr(trms, "term.labels")
    M = model.frame(trms, utils::head(data, 2)) #### just need a couple rows
    isfac = sapply(M, function(x) inherits(x, "factor"))
    
    # Character vector of terms in the model frame that are factors ...
    facs.m = names(M)[isfac]
    
    # Exclude the terms that are already factors
    # What's left will be things like "factor(dose)", "interact(dose,treat)", etc
    cterms = setdiff(facs.m, facs.d)
    
    if(length(cterms) == 0) 
        return(cterms)
    # (else) Strip off the function calls
    cvars = lapply(cterms, function(x) .all.vars(reformulate(x)))
    
    # Exclude any variables that are already factors
    intersect(unique(unlist(cvars)), covs.d)
}

# calculate the offset for the given grid
.get.offset = function(terms, grid) {
    off.idx = attr(terms, "offset")
    offset = rep(0, nrow(grid))
    tvars = attr(terms, "variables")
    for (i in off.idx)
        offset = offset + eval(tvars[[i+1]], grid)
    offset
}



### =========== Methods for ref.grid class =============================
# (note: summary-related methods moved to a new file)

str.ref.grid <- function(object, ...) {
    showlevs = function(x) { # internal convenience function
        if (is.null(x)) cat("(predicted by other variables)")
        else cat(paste(format(x, digits = 5, justify = "none"), collapse=", "))
    }
    #cat("responses: ")
    #showlevs(object@roles$responses)
    levs = object@levels
    cat(paste("'", class(object)[1], "' object with variables:\n", sep=""))
    for (nm in union(object@roles$predictors, union(object@roles$multresp, object@roles$responses))) {
        cat(paste("    ", nm, " = ", sep = ""))
        if (nm %in% names(object@matlevs)) {
            if (nm %in% object@roles$responses)
                cat("multivariate response with means: ")
            else
                cat("matrix with column means: ")
            cat("\n        ")
            showlevs(object@matlevs[[nm]])
        }
        else if (nm %in% object@roles$multresp) {
            cat("multivariate response levels: ")
            showlevs(levs[[nm]])
        }
        else if (nm %in% object@roles$responses) {
            cat("response variable with mean ")
            showlevs(levs[[nm]])
        }
        else
            showlevs(levs[[nm]])
        cat("\n")
    }
    if(!is.null(tran <- object@misc$tran)) {
        if (is.list(tran)) 
            tran = ifelse(is.null(tran$name), "custom - see slot(, \"misc\")$tran", tran$name)
        if (!is.null(mul <- object@misc$tran.mult))
            tran = paste0(mul, "*", tran)
        cat(paste("Transformation:", dQuote(tran), "\n"))
    }
}



print.ref.grid = function(x,...)
    print(summary.ref.grid(x, ...))


# vcov method
vcov.ref.grid = function(object, ...) {
    tol = get.lsm.option("estble.tol")
    if (!is.null(hook <- object@misc$vcovHook)) {
        if (is.character(hook)) 
            hook = get(hook)
        hook(object, tol = tol, ...)
    }
    else {
        X = object@linfct
        estble = estimability::is.estble(X, object@nbasis, tol) ###apply(X, 1, .is.estble, object@nbasis, tol)
        X[!estble, ] = NA
        X = X[, !is.na(object@bhat), drop = FALSE]
        X %*% tcrossprod(object@V, X)
    }
}


# Method to alter contents of misc slot
update.ref.grid = function(object, ..., silent = FALSE) {
    args = list(...)
    valid.misc = c("adjust","alpha","avgd.over","by.vars","delta","df",
        "initMesg","estName","estType","famSize","infer","inv.lbl",
        "level","methdesc","null","predict.type","pri.vars","side","tran","tran.mult")
    valid.slots = slotNames(object)
    valid.choices = union(valid.misc, valid.slots)
    misc = object@misc
    for (nm in names(args)) {
        fullname = try(match.arg(nm, valid.choices), silent=TRUE)
        if(inherits(fullname, "try-error")) {
            if (!silent)
                message("Argument ", sQuote(nm), " was ignored. Valid choices are:\n",
                    paste(valid.choices, collapse=", "))
        }
        else {
            if (fullname %in% valid.slots)
                slot(object, fullname) = args[[nm]]
            else {
                if (fullname == "by.vars") {
                    allvars = union(misc$pri.vars, misc$by.vars)
                    misc$pri.vars = setdiff(allvars, args[[nm]])
                }
                if (fullname == "pri.vars") {
                    allvars = union(misc$pri.vars, misc$by.vars)
                    misc$by.vars = setdiff(allvars, args[[nm]])
                }
                misc[[fullname]] = args[[nm]]
            }
        }
    }
    object@misc = misc
    object
}

### set or change lsmeans options
lsm.options = function(...) {
    opts = getOption("lsmeans", list())
#    if (is.null(opts)) opts = list()
    newopts = list(...)
    for (nm in names(newopts))
        opts[[nm]] = newopts[[nm]]
    options(lsmeans = opts)
    invisible(opts)
}

# equivalent of getOption()
get.lsm.option = function(x, default = defaults.lsm[[x]]) {
    opts = getOption("lsmeans", list())
    if(is.null(default) || x %in% names(opts))
        opts[[x]]
    else 
        default
}

### Exported defaults for certain options
defaults.lsm = list(
    estble.tol = 1e-8,        # tolerance for estimability checks
    disable.pbkrtest = FALSE, # whether to bypass pbkrtest routines for lmerMod
    pbkrtest.limit = 3000,    # limit on N for enabling adj V
    save.ref.grid = TRUE      # save new ref.grid in .Last.ref.grid
)

# Utility that returns TRUE if getOption("lsmeans")[[opt]] is TRUE
.lsm.is.true = function(opt) {
    x = get.lsm.option(opt, FALSE)
    if (is.logical(x))  x
    else FALSE
}


### Utility to change the internal structure of a ref.grid
### Returned ref.grid object has linfct = I and bhat = estimates
### Primary reason to do this is with transform = TRUE, then can 
### work with linear functions of the transformed predictions
regrid = function(object, transform = c("response", "log", "none"), inv.log.lbl = "response") {
    if (is.logical(transform))   # for backward-compatibility
        transform = ifelse(transform, "response", "none")
    else
        transform = match.arg(transform)
    est = .est.se.df(object, do.se = TRUE) ###FALSE)
    estble = !(is.na(est[[1]]))
    object@V = vcov(object)[estble, estble, drop=FALSE]
    object@bhat = est[[1]]
    object@linfct = diag(1, length(estble))
    if(all(estble))
        object@nbasis = estimability::all.estble
    else
        object@nbasis = object@linfct[, !estble, drop = FALSE]
    
    # override the df function
    df = est$df
    test.df = diff(range(df))
    if (is.na(test.df) || test.df < .001) {
        object@dfargs = list(df = mean(df))
        object@dffun = function(k, dfargs) dfargs$df
    }
    else { # use containment df
        object@dfargs = list(df = df)
        object@dffun = function(k, dfargs) {
            idx = which(zapsmall(k) != 0)
            ifelse(length(idx) == 0, NA, min(dfargs$df[idx]))
        }
    }
    
    if(transform %in% c("response", "log") && !is.null(object@misc$tran)) {
        link = attr(est, "link")
        D = .diag(link$mu.eta(object@bhat[estble]))
        object@bhat = sapply(object@bhat, function(x) 
            ifelse(link$valideta(x), link$linkinv(x), 0))
        object@V = D %*% tcrossprod(object@V, D)
        inm = object@misc$inv.lbl
        if (!is.null(inm))
            object@misc$estName = inm
        object@misc$tran = object@misc$tran.mult = object@misc$inv.lbl = NULL
    }
    if (transform == "log") {
        D = .diag(1/object@bhat)
        object@V = D %*% tcrossprod(object@V, D)
        object@bhat = log(object@bhat)
        object@misc$tran = "log"
        object@misc$inv.lbl = inv.log.lbl
    }
    # Nix out things that are no longer needed or valid
    object@grid$.offset. = object@misc$offset.mult =
        object@misc$estHook = object@misc$vcovHook = NULL
    object
}


### S4 methods
## use S3 for this setMethod("summary", "ref.grid", summary.ref.grid)
setMethod("show", "ref.grid", function(object) str.ref.grid(object))



