### Code for an enhancement of 'glht' in 'multcomp' package
### Provides for using 'lsm' in similar way to 'mcp'
### This is implemented via the class "lsmlf" -- linear functions for lsmeans
### (also oddly reminiscent of an old Lucky Strike commercial, LSMFT)

# lsm(specs) will be used as 'linfct' argument in glht
# all we need to do is class it and save the arguments
lsm <- function(...) {
    result <- list(...)
    class(result) <- "lsmlf"
    result
}

# New S3 method for lsmlf objects
glht.lsmlf <- function(model, linfct, ...) {
    # Pass the arguments we should pass to ref.grid:
    args = linfct
    args[[1]] = model
    names(args)[1] = "object"
    # Now pass the ref.grid to lsmeans:
    linfct$object <- do.call("ref.grid", args)
    lsmo <- do.call("lsmeans", linfct)
    if (is.list(lsmo)) 
        lsmo = lsmo[[length(lsmo)]]
    # Then call the method for lsmobj
    glht(model, lsmo, ...)
}


# S3 method for an lsmobj or ref.grid
# Note: model is redundant, really, so can be omitted
glht.ref.grid <- function(model, linfct, by, ...) {
    if (!requireNamespace("multcomp"))
        stop(sQuote("glht")," requires ", dQuote("multcomp"), " to be installed")
    object = linfct # so I don't get confused
    if (missing(model)) 
        model = .cls.list("lsmwrap", object = object)
    args = list(model = model, ...)
    # add a df value if not supplied
    if (is.null(args$df)) {
        df = summary(linfct)$df
        if(any(!is.na(df))) {
            args$df = max(1, as.integer(mean(df, na.rm=TRUE) + .25))
            message("Note: df set to ", args$df)
        }
    }
    if (missing(by)) by = object@misc$by.vars
    
    nms = setdiff(names(object@grid), c(by, ".offset.", ".freq.", ".wgt."))
    if (is.null(object@misc$estHook))
        lf = object@linfct
    else # custom estimation setup - use the grid itself as the parameterization
        lf = diag(1, nrow(object@linfct))
    dimnames(lf)[[1]] = as.character(interaction(object@grid[, nms], sep=", "))
    
    if (is.null(by)) {
        args$linfct = lf
        return(do.call("glht", args))
    }
    
    # (else...)
    by.rows = .find.by.rows(object@grid, by)
    result = lapply(by.rows, function(r) {
        args$linfct = lf[r, , drop=FALSE]
        do.call("glht", args)
    })
    bylevs = lapply(by, function(byv) unique(object@grid[[byv]]))
    names(bylevs) = by
    bygrid = do.call("expand.grid", bylevs)
    levlbls = lapply(by, function(byv) paste(byv, "=", bygrid[[byv]]))
    levlbls$sep = ", "
    names(result) = do.call("paste", levlbls)
    class(result) = c("glht.list", "list")
    result
}

### as. glht -- convert my object to glht object
as.glht <- function(object, ...)
    UseMethod("as.glht")

as.glht.default <- function(object, ...)
    stop("Cannot convert an object of class ", sQuote(class(object)[1]),
         " to a ", sQuote("glht"), " object")

as.glht.ref.grid <- function(object, ...)
    glht( , object, ...)

as.glht.lsm.list <- function(object, ..., which = 1)
    as.glht(object[[which]], ...)


# S3 modelparm method for lsmwrap (S3 wrapper for an lsmobj - see glht.lsmobj)
modelparm.lsmwrap <- function(model, coef., vcov., df, ...) {
    object = model$object
    if (is.null(object@misc$estHook)) {
        bhat = object@bhat
        V = object@V
    }
    else { # Have custom vcov and est methods. Use the grid itself as parameterization
        bhat = predict(object)
        V = vcov(object)
    }
    if(missing(df) || is.na(df))
        df = 0
    .cls.list("modelparm", coef = bhat, vcov = V,
                df = df, estimable = !is.na(bhat))
    # This is NOT what we mean by 'estimable', but it is what glht wants...
}

# S3 methods for glht.list

### Doesn't work so excluded...
# cld.glht.list = function(object, ...)
#     lapply(object, cld, ...)

coef.glht.list = function(object, ...)
    lapply(object, coef, ...)

confint.glht.list = function(object, ...)
    lapply(object, confint, ...)

plot.glht.list = function(x, ...)
    lapply(x, plot, ...)

summary.glht.list = function(object, ...)
    lapply(object, summary, ...)

vcov.glht.list = function(object, ...)
    lapply(object, vcov, ...)






