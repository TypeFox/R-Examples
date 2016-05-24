
# $Id: glht.R 422 2015-07-20 13:13:04Z thothorn $

### general linear hypotheses
glht <- function(model, linfct, ...) {
    if (missing(linfct)) {
        mpar <- modelparm(model, ...)
        linfct <- diag(length(mpar$coef))
        rownames(linfct) <- names(mpar$coef)
        glht(model = model, linfct = linfct, ...)
    } else {
        UseMethod("glht", linfct)
    }
}

### K coef(model) _!alternative_ rhs
glht.matrix <- function(model, linfct, 
    alternative = c("two.sided", "less", "greater"), rhs = 0, ...) {

    ### extract coefficients and their covariance matrix, df
    mpar <- modelparm(model, ...)

    alternative <- match.arg(alternative)
    if (!is.numeric(rhs))
        stop(sQuote("rhs"), " is not a numeric vector")

    if (!all(mpar$estimable)) {
        ignoreOK <- all(colSums(abs(linfct[, !mpar$estimable, 
                                           drop = FALSE])) == 0)
        if (!ignoreOK)
            stop("some linear functions are not estimable")
        linfct <- linfct[, mpar$estimable, drop = FALSE]
        warning(sum(!mpar$estimable), " out of ", length(mpar$estimable), 
                " coefficients not estimable in ", sQuote("model"))
    }

    if (ncol(linfct) != length(mpar$coef))
        stop(sQuote("ncol(linfct)"), " is not equal to ", 
             sQuote("length(coef(model))"))

    if (is.null(colnames(linfct)))
        colnames(linfct) <- names(mpar$coef)

    if (is.null(rownames(linfct))) # {
        rownames(linfct) <- 1:nrow(linfct)
#    } else {
        ### alt <- switch(alternative, 
        ###    "two.sided" = "==", "less" = ">=", "greater" = "<=")
        ### rownames(linfct) <- paste(rownames(linfct), alt, rhs)
#    }

    if (length(rhs) == 1) rhs <- rep(rhs, nrow(linfct))
    if (length(rhs) != nrow(linfct))
        stop(sQuote("nrow(linfct)"), " is not equal to ",
             sQuote("length(rhs)"))

    RET <- list(model = model, linfct = linfct, rhs = rhs,
                coef = mpar$coef, vcov = mpar$vcov, 
                df = mpar$df, alternative = alternative,
                type = attr(linfct, "type"))
    class(RET) <- "glht"
    RET
}

### symbolic description of H_0
glht.character <- function(model, linfct, ...) {

    ### extract coefficients and their covariance matrix
    beta <- modelparm(model, ...)$coef
    tmp <- chrlinfct2matrix(linfct, names(beta))
    return(glht(model, linfct = tmp$K, rhs = tmp$m, 
                alternative = tmp$alternative, ...))
}

### symbolic description of H_0
glht.expression <- function(model, linfct, ...) 
    glht(model, deparse(linfct), ...)

### multiple comparison procedures
glht.mcp <- function(model, linfct, ...) {

    ### extract factors and contrast matrices from `model'
    ia <- attr(linfct, "interaction_average")
    ca <- attr(linfct, "covariate_average")
    if (ia || ca) {
        ### experimental version
        tmp <- mcp2matrix2(model, linfct = linfct, interaction_average = ia,
                           covariate_average = ca)
    } else {
        ### use old version
        tmp <- mcp2matrix(model, linfct = linfct)
    }
    args <- list(model = model, linfct = tmp$K)
    if (!is.null(tmp$alternative))
        args$alternative <- tmp$alternative
    if (any(tmp$m != 0))
        args$rhs <- tmp$m
    args <- c(args, list(...))

    ret <- do.call("glht", args)
    ret$type <- tmp$type
    ret$focus <- names(linfct)
    return(ret)
}

### call Rich' function for raw means ...
glht.means <- function(model, linfct, ...) {

    args <- list(model = model, 
                 linfct = meanslinfct(model, focus = names(linfct), ...))
    args <- c(args, list(...))
    ret <- do.call("glht", args)
    ret$type <- "Mean"
    ret$focus <- names(linfct)
    return(ret)
}
