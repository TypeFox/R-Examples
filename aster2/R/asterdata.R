
asterdata <- function(data, vars, pred, group, code, families, delta,
    response.name = "resp", varb.name = "varb",
    tolerance = 8 * .Machine$double.eps) {
    stopifnot(is.data.frame(data))
    stopifnot(is.atomic(vars))
    stopifnot(is.character(vars))
    stopifnot(all(vars %in% names(data)))
    if ("id" %in% names(data))
        warning("The reshape function (which this function uses) by default",
            "\ntakes a column named \"id\" as its idvar.  Is this wanted?")
    stopifnot(is.atomic(pred))
    stopifnot(is.numeric(pred))
    stopifnot(length(pred) == length(vars))
    stopifnot(all(pred == as.integer(pred)))
    stopifnot(all(pred >= 0))
    stopifnot(all(pred < seq(along = pred)))
    stopifnot(is.atomic(group))
    stopifnot(is.numeric(group))
    stopifnot(length(group) == length(vars))
    stopifnot(all(group == as.integer(group)))
    stopifnot(all(group >= 0))
    stopifnot(all(group < seq(along = group)))
    pred.group <- c(0, pred)[group + 1]
    if (! all(group == 0 | pred.group == pred))
       stop("variables in some dependence group have different predecessors")
    stopifnot(is.atomic(code))
    stopifnot(is.numeric(code))
    stopifnot(length(code) == length(vars))
    stopifnot(all(code == as.integer(code)))
    stopifnot(all(code %in% seq(along = families)))
    stopifnot(is.list(families))
    if (! all(sapply(families, function(x) inherits(x, "astfam") |
        (is.character(x) & length(x) == 1))))
        stop("some component of families not class \"astfam\" or string")
    for (i in seq(along = families))
        if (is.character(families[[i]])) {
            foo <- try(eval(call(paste("fam", families[[i]], sep = "."))))
            if (inherits(foo, "try-error"))
                stop("invalid character string abbreviation for family")
            families[[i]] <- foo
        }
    stopifnot(all(sapply(families, function(x) inherits(x, "astfam"))))
    if (missing(delta)) delta <- rep(0, length(vars))
    stopifnot(is.atomic(delta))
    stopifnot(is.numeric(delta))
    stopifnot(length(delta) == length(vars))
    stopifnot(all(is.finite(delta)))
    stopifnot(is.atomic(response.name))
    stopifnot(is.character(response.name))
    stopifnot(length(response.name) == 1)
    stopifnot(is.atomic(varb.name))
    stopifnot(is.character(varb.name))
    stopifnot(length(varb.name) == 1)
    stopifnot(is.atomic(tolerance))
    stopifnot(is.numeric(tolerance))
    stopifnot(length(tolerance) == 1)
    stopifnot(tolerance > 0)
    ### at this point arguments are valid except we still need to check that
    ### each dependence group has the correct dimension
    fam.clear()
    for (i in seq(along = families)) fam.set(families[[i]])
    d <- sapply(as.list(code), fam.dimension)
    group.idx <- seq(along = vars)
    repeat {
        group.idx.save <- group.idx
        group.idx[group > 0] <- group.idx[group]
        if (identical(group.idx, group.idx.save)) break
    }
    group.dimen <- tabulate(group.idx, nbins = length(code))
    if (! all(group.dimen == 0 | group.dimen == d)) 
       stop("number of variables in some dependence group not equal dimension")
    ### at this point arguments are valid and families have been set up
    redata <- reshape(data, varying = list(vars), direction = "long",
        timevar = varb.name, times = factor(vars, levels = vars),
        v.names = response.name)
    repred <- rep(pred, each = nrow(data))
    repred <- (repred - 1) * nrow(data)
    repred <- repred + rep(1:nrow(data), times = length(vars))
    repred[repred < 0] <- 0
    regroup <- rep(group, each = nrow(data))
    regroup <- (regroup - 1) * nrow(data)
    regroup <- regroup + rep(1:nrow(data), times = length(vars))
    regroup[regroup < 0] <- 0
    recode <- rep(code, each = nrow(data))
    redelta <- rep(delta, each = nrow(data))
    fam.clear()
    result <- structure(list(redata = redata, repred = repred,
        regroup = regroup, recode = recode, families = families,
        redelta = redelta, initial = rep(1, nrow(redata)),
        response.name = response.name, pred = pred, group = group,
        code = code), class = "asterdata")
    validasterdata(result, tolerance)
    return(result)
}

validasterdata <- function(object, tolerance = 8 * .Machine$double.eps) {
    stopifnot(inherits(object, "asterdata"))
    stopifnot(is.list(object))
    stopifnot(all(c("redata", "repred", "initial", "regroup", "recode",
        "families", "redelta", "response.name") %in% names(object)))
    stopifnot(is.atomic(tolerance))
    stopifnot(is.numeric(tolerance))
    stopifnot(length(tolerance) == 1)
    stopifnot(tolerance > 0)
    # response.name
    response.name <- object$response.name
    stopifnot(is.atomic(response.name))
    stopifnot(is.character(response.name))
    stopifnot(length(response.name) == 1)
    # redata and resp
    stopifnot(is.data.frame(object$redata))
    stopifnot(response.name %in% names(object$redata))
    resp <- object$redata[[response.name]]
    # repred
    repred <- object$repred
    stopifnot(is.atomic(repred))
    stopifnot(is.numeric(repred))
    stopifnot(length(repred) == length(resp))
    stopifnot(all(repred == as.integer(repred)))
    stopifnot(all(repred >= 0))
    stopifnot(all(repred < seq(along = repred)))
    # regroup
    regroup <- object$regroup
    stopifnot(is.atomic(regroup))
    stopifnot(is.numeric(regroup))
    stopifnot(length(regroup) == length(resp))
    stopifnot(all(regroup == as.integer(regroup)))
    stopifnot(all(regroup >= 0))
    stopifnot(all(regroup < seq(along = regroup)))
    # families
    families <- object$families
    stopifnot(is.list(families))
    stopifnot(all(sapply(families, function(x) inherits(x, "astfam"))))
    # recode
    recode <- object$recode
    stopifnot(is.atomic(recode))
    stopifnot(is.numeric(recode))
    stopifnot(length(recode) == length(resp))
    stopifnot(all(recode == as.integer(recode)))
    stopifnot(all(recode %in% seq(along = families)))
    # redelta
    redelta <- object$redelta
    stopifnot(is.atomic(redelta))
    stopifnot(is.numeric(redelta))
    stopifnot(length(redelta) == length(resp))
    stopifnot(all(is.finite(redelta)))
    # initial
    initial <- object$initial
    stopifnot(is.atomic(initial))
    stopifnot(is.numeric(initial))
    stopifnot(length(initial) == length(resp))
    stopifnot(all(is.finite(initial)))
    stopifnot(all(initial > 0))
    ### at this point all components are valid except we still need to check
    ### that
    ### (1) each dependence group has the correct dimension
    ### (2) predecessor values are valid
    ### (3) delta is valid
    ### (4) response values are valid, given predecessor and delta
    fam.set.tolerance(tolerance)
    # set up families
    fam.clear()
    for (i in seq(along = families)) fam.set(families[[i]])
    .C("aster_validate", nnode = length(resp), resp = as.double(resp),
       pred = as.integer(repred), group = as.integer(regroup),
       code = as.integer(recode), initial = as.double(initial),
       delta = as.double(redelta), PACKAGE = "aster2")
    fam.clear()
    fam.reset.tolerance()
    invisible(TRUE)
}

is.validasterdata <- function(object, tolerance = 8 * .Machine$double.eps) {
    stopifnot(inherits(object, "asterdata"))
    old.opt <- options(show.error.messages = FALSE)
    out <- try(validasterdata(object, tolerance))
    options(old.opt)
    return(! inherits(out, "try-error"))
}

subset.asterdata <- function(x, subset, successors = TRUE, ...) {
    stopifnot(inherits(x, "asterdata"))
    validasterdata(x)
    stopifnot(is.logical(subset))
    stopifnot(length(subset) == nrow(x$redata))
    stopifnot(is.logical(successors))
    stopifnot(length(successors) == 1)

    idx <- seq(along = subset)[subset]

    regroup.sub.gpred <- x$regroup[subset]
    regroup.sub.gpred <- regroup.sub.gpred[regroup.sub.gpred != 0]
    regroup.sub.gsucc <- seq(along = x$regroup)[x$regroup %in% idx]
    if ((! all(regroup.sub.gpred %in% idx)) ||
        (! all(regroup.sub.gsucc %in% idx)))
        stop("some dependence group not all in or all out of subset")

    repred.sub.gpred <- x$repred[subset]
    repred.sub.gpred <- repred.sub.gpred[repred.sub.gpred != 0]
    repred.sub.gsucc <- seq(along = x$repred)[x$repred %in% idx]
    if ((! all(repred.sub.gpred %in% idx)))
        stop("some predecessors of subset not in subset")
    if (successors && (! all(repred.sub.gsucc %in% idx)))
        stop("some successors of subset not in subset")

    subset[is.na(subset)] <- FALSE

    new.idx <- seq(along = idx)
    new.redata <- subset(x$redata, subset)
    new.repred <- x$repred[subset]
    new.repred <- match(new.repred, idx)
    new.repred[is.na(new.repred)] <- 0
    new.regroup <- x$regroup[subset]
    new.regroup <- match(new.regroup, idx)
    new.regroup[is.na(new.regroup)] <- 0
    result <- list(redata = new.redata, repred = new.repred,
        regroup = new.regroup, recode = x$recode[subset],
        families = x$families, redelta = x$redelta[subset],
        initial = x$initial[subset], response.name = x$response.name)
    result$pred <- x$pred
    result$group <- x$group
    result$code <- x$code
    class(result) <- "asterdata"
    validasterdata(result)
    return(result)
}

