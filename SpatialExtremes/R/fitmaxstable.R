fitmaxstab <- function(data, coord, cov.mod, loc.form, scale.form, shape.form,
                       marg.cov = NULL, temp.cov = NULL, temp.form.loc = NULL,
                       temp.form.scale = NULL, temp.form.shape = NULL, iso = FALSE,
                       ..., fit.marge = FALSE, warn = TRUE, method = "Nelder", start,
                       control = list(), weights = NULL, corr = FALSE, check.grad = FALSE){

    if (is.null(dim(coord))){
        if (length(coord) != ncol(data))
            stop("'data' and 'coord' don't match")
    }

    else if (nrow(coord) != ncol(data))
        stop("'data' and 'coord' don't match")

    if (!is.null(marg.cov) && is.null(colnames(marg.cov)))
        stop("'marg.cov' must have named columns")

    if (!is.null(marg.cov) && (nrow(marg.cov) != nrow(coord)))
        stop("'data' and 'marg.cov' don't match")

    if (missing(loc.form) && missing(scale.form) && missing(shape.form))
        reg.mod <- "full"

    if (!missing(loc.form) && !missing(scale.form) && !missing(shape.form)){
        reg.mod <- "spatgev"
        fit.marge <- TRUE

        if ((class(loc.form) != "formula") || (class(scale.form) != "formula") ||
            (class(shape.form) != "formula"))
            stop("''loc.form'', ''scale.form'' and ''shape.form'' must be valid R formulas")
    }

    flag <- missing(loc.form) + missing(scale.form)  + missing(shape.form)

    if (!(flag %in% c(0, 3)))
        stop("if one formula is given for the GEV parameters, then it should
be given for *ALL* GEV parameters")

    n.pairs <- ncol(data) * (ncol(data) - 1) / 2
    if (!is.null(weights) && (!is.numeric(weights) || length(weights) != n.pairs ||
                              all(weights == 0)))
        stop("The weights you specified are not valid")

    if (method != "nlminb"){
        if (is.null(control$maxit))
            control$maxit <- 10000
    }

    else{
        if (is.null(control$eval.max))
            control$eval.max <- 15000

        if (is.null(control$iter.max))
            control$iter.max <- 10000
    }

    if (check.grad)
        ## Force to use the nlm optimizer as it returns the numerical gradient
        method <- "nlm"

    if (cov.mod == "gauss")
        fitted <- switch(reg.mod,
                         "full" = smithfull(data, coord, ..., fit.marge = fit.marge, iso = iso,
                         warn = warn, method = method, control = control, check.grad = check.grad,
                         corr = corr, start = start, weights = weights),
                         "spatgev" = smithform(data, coord, ..., loc.form = loc.form, scale.form = scale.form,
                         shape.form = shape.form, fit.marge = fit.marge, iso = iso, marg.cov = marg.cov,
                         warn = warn, method = method, control = control, corr = corr, start = start, weights = weights,
                         temp.cov = temp.cov, temp.form.loc = temp.form.loc, check.grad = check.grad,
                         temp.form.scale = temp.form.scale, temp.form.shape = temp.form.shape))

    else if (cov.mod == "brown")
        fitted <- switch(reg.mod,
                         "full" = brownresnickfull(data, coord, ..., fit.marge = fit.marge,
                         warn = warn, method = method, control = control,
                         corr = corr, start = start, weights = weights, check.grad = check.grad),
                         "spatgev" = brownresnickform(data, coord, ..., loc.form = loc.form,
                         scale.form = scale.form, shape.form = shape.form, check.grad = check.grad,
                         fit.marge = fit.marge, marg.cov = marg.cov, warn = warn,
                         method = method, control = control, corr = corr,
                         start = start, weights = weights, temp.cov = temp.cov, temp.form.loc = temp.form.loc,
                         temp.form.scale = temp.form.scale, temp.form.shape = temp.form.shape))
    else if (substr(cov.mod, 1, 1) == "i")
        fitted <- switch(reg.mod,
                         "full" = schlatherindfull(data, coord, cov.mod = substr(cov.mod, 2, 8),
                         ..., fit.marge = fit.marge, warn = warn, check.grad = check.grad,
                         method = method, control = control,
                         corr = corr, start = start, weights = weights),
                         "spatgev" = schlatherindform(data, coord, cov.mod = substr(cov.mod, 2, 8), ...,
                         loc.form = loc.form, scale.form = scale.form, shape.form = shape.form,
                         fit.marge = fit.marge, marg.cov = marg.cov, warn = warn, check.grad = check.grad,
                         method = method, control = control, corr = corr,
                         start = start, weights = weights, temp.cov = temp.cov, temp.form.loc = temp.form.loc,
                         temp.form.scale = temp.form.scale, temp.form.shape = temp.form.shape))

    else if (substr(cov.mod, 1, 1) == "g")
        fitted <- switch(reg.mod,
                         "full" = geomgaussfull(data, coord, cov.mod = substr(cov.mod, 2, 8),
                         ..., fit.marge = fit.marge, warn = warn, check.grad = check.grad,
                         method = method, control = control,
                         corr = corr, start = start, weights = weights),
                         "spatgev" = geomgaussform(data, coord, cov.mod = substr(cov.mod, 2, 8), ...,
                         loc.form = loc.form, scale.form = scale.form, shape.form = shape.form,
                         fit.marge = fit.marge, marg.cov = marg.cov, warn = warn, check.grad = check.grad,
                         method = method, control = control, corr = corr,
                         start = start, weights = weights, temp.cov = temp.cov, temp.form.loc = temp.form.loc,
                         temp.form.scale = temp.form.scale, temp.form.shape = temp.form.shape))

    else if (substr(cov.mod, 1, 1) == "t")
        fitted <- switch(reg.mod,
                         "full" = extremaltfull(data, coord, cov.mod = substr(cov.mod, 2, 8),
                         ..., fit.marge = fit.marge, warn = warn, check.grad = check.grad,
                         method = method, control = control,
                         corr = corr, start = start, weights = weights),
                         "spatgev" = extremaltform(data, coord, cov.mod = substr(cov.mod, 2, 8), ...,
                         loc.form = loc.form, scale.form = scale.form, shape.form = shape.form,
                         fit.marge = fit.marge, marg.cov = marg.cov, warn = warn, check.grad = check.grad,
                         method = method, control = control, corr = corr,
                         start = start, weights = weights, temp.cov = temp.cov, temp.form.loc = temp.form.loc,
                         temp.form.scale = temp.form.scale, temp.form.shape = temp.form.shape))

    else
        fitted <- switch(reg.mod,
                         "full" = schlatherfull(data, coord, cov.mod = cov.mod,
                         ..., fit.marge = fit.marge, warn = warn, check.grad = check.grad,
                         method = method, control = control,
                         corr = corr, start = start, weights = weights),
                         "spatgev" = schlatherform(data, coord, cov.mod = cov.mod, ...,
                         loc.form = loc.form, scale.form = scale.form, shape.form = shape.form,
                         fit.marge = fit.marge, marg.cov = marg.cov, warn = warn, check.grad = check.grad,
                         method = method, control = control, corr = corr,
                         start = start, weights = weights, temp.cov = temp.cov, temp.form.loc = temp.form.loc,
                         temp.form.scale = temp.form.scale, temp.form.shape = temp.form.shape))

    return(fitted)
}

##fitnsmaxstab## <- function(data, coord, cov.mod, sigma2.form, loc.form, scale.form, shape.form,
##                         marg.cov = NULL, ..., fit.marge = FALSE, warn = TRUE,
##                         method = "Nelder", start, control = list(),
##                         std.err.type = "score", corr = FALSE){
##
##  if (!(std.err.type) %in% c("none", "score", "grad"))
##    stop("'std.err.type' must be one of 'none', 'score' or 'grad'")
##
##  if (nrow(coord) != ncol(data))
##    stop("'data' and 'coord' don't match")
##
##  if (!is.null(marg.cov) && is.null(colnames(marg.cov)))
##    stop("'marg.cov' must have named columns")
##
##  if (!is.null(marg.cov) && (nrow(marg.cov) != nrow(coord)))
##    stop("'data' and 'marg.cov' don't match")
##
##  if (missing(loc.form) && missing(scale.form) && missing(shape.form))
##    reg.mod <- "full"
##
##  if (!missing(loc.form) && !missing(scale.form) && !missing(shape.form)){
##    reg.mod <- "spatgev"
##    fit.marge <- TRUE
##
##    if ((class(loc.form) != "formula") || (class(scale.form) != "formula") ||
##        (class(shape.form) != "formula"))
##      stop("''loc.form'', ''scale.form'' and ''shape.form'' must be valid R formulas")
##  }
##
##  flag <- missing(loc.form) + missing(scale.form)  + missing(shape.form)
##
##  if (!(flag %in% c(0, 3)))
##    stop("if one formula is given for the GEV parameters, then it should
##be given for *ALL* GEV parameters")
##
##  if (method != "nlminb"){
##    if (is.null(control$maxit))
##      control$maxit <- 10000
##  }
##
##  else{
##    if (is.null(control$eval.max))
##      control$eval.max <- 15000
##
##    if (is.null(control$iter.max))
##      control$iter.max <- 10000
##  }
##
##  fitted <- switch(reg.mod, "full" = nsgeomgaussfull(data, coord, cov.mod = cov.mod,
##                              sigma2.form, ..., fit.marge = fit.marge, warn = warn,
##                              method = method, control = control, std.err.type = std.err.type,
##                              corr = corr, start = start),
##                   "spatgev" = nsgeomgaussform(data, coord, cov.mod = cov.mod,
##                     sigma2.form, ..., loc.form = loc.form, scale.form = scale.form,
##                     shape.form = shape.form, fit.marge = fit.marge, marg.cov = marg.cov, warn = warn,
##                     method = method, control = control, std.err.type = std.err.type, corr = corr,
##                     start = start))
##
##  return(fitted)
##}
