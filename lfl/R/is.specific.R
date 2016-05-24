.is.specific <- function(x, y, vars, specs) {
    lvl <- names(vars)

    config <- list(x=as.integer(factor(x, levels=lvl)) - 1,
                   y=as.integer(factor(y, levels=lvl)) - 1,
                   vars=as.numeric(as.factor(vars[lvl])),
                   specs=specs[lvl, lvl])
    result <- .Call("isSpecific", config, PACKAGE="lfl")
    return(result)
}


# Determines whether 'x' is more specific (or equal!!!) than 'y' with respect to 'vars' and 'specs'.
is.specific <- function(x, y, vars, specs) {
    if (!is.null(x) && (!is.vector(x) || !is.character(x))) {
        stop("'x' must to be a character vector")
    }
    if (!is.null(y) && (!is.vector(y) || !is.character(y))) {
        stop("'y' must to be a character vector")
    }
    if (!is.vector(vars) || is.null(names(vars))) {
        stop("'vars' must be a named vector")
    }

    xVars <- vars[x]
    yVars <- vars[y]

    if (any(is.na(c(xVars, yVars)))) {
        stop("'vars' not compatible with input 'x' or 'y'")
    }
    if (!all(sort(names(vars)) == sort(colnames(specs))) ||
            !all(sort(names(vars)) == sort(colnames(specs)))) {
        stop("'specs' colnames or rownames are incompatible with 'vars' names")
    }
    if ((length(unique(xVars)) != length(xVars)) || (length(unique(yVars)) != length(yVars))) {
        stop('Unable to work with rules containing the same var more times')
    }

    return(.is.specific(x, y, vars, specs))
}
