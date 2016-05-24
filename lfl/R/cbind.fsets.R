cbind.fsets <- function(..., deparse.level = 1) {
    dots <- list(...)

    res <- NULL
    resVars <- NULL
    resSpecs <- NULL

    for (i in seq_along(dots)) {
        arg <- dots[[i]]
        argName <- names(dots)[i]

        if (!is.null(arg)) {
            if (!is.fsets(arg)) {
                stop("Function 'cbind.fsets' cannot bind arguments that are not valid 'fsets' objects")
            }
            class(arg) <- setdiff(class(arg), 'fsets')
            if (is.null(res)) {
                resVars <- vars(arg)
                resSpecs <- specs(arg)
                res <- arg
            } else {
                resVarNames <- c(names(resVars), names(vars(arg)))
                resVars <- c(resVars, vars(arg))
                names(resVars) <- resVarNames

                o1 <- matrix(0, nrow=nrow(resSpecs), ncol=ncol(specs(arg)))
                o2 <- matrix(0, nrow=nrow(specs(arg)), ncol=ncol(resSpecs))
                resSpecs <- rbind(cbind(resSpecs, o1),
                                  cbind(o2, specs(arg)))
                colnames(resSpecs) <- names(resVars)
                rownames(resSpecs) <- names(resVars)
                res <- cbind(res, arg)
            }
        }
    }

    return(fsets(res, resVars, resSpecs))
}
