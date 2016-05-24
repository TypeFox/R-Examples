
as.basis.factor_var <- function(object, ...) {
    fm <- as.formula(paste("~ 0 +", variable.names(object)))
    nd <- data.frame(support(object))
    names(nd) <- variable.names(object)
    as.basis(fm, data = nd, ...)
}

as.basis.ordered_var <- function(object, ...) {
    fm <- as.formula(paste("~ ", variable.names(object)))
    nl <- nlevels(object)
    stopifnot(nl > 2)
    nd <- data.frame(support(object))
    names(nd) <- variable.names(object)
    ctr <- list(function(n) contr.treatment(n, base = nl))
    names(ctr) <- variable.names(object)
    as.basis(fm, data = nd, remove_intercept = TRUE,
             contrasts.arg = ctr, ui = diff(diag(nl - 1)), 
             ci = rep(0, nl - 2))
} 

as.basis.factor <- function(object, ...) {
    vn <- deparse(substitute(object))
    vn <- gsub(".*\\$", "", vn)
    as.basis(factor_var(vn, levels = levels(object)))
}

as.basis.ordered <- function(object, ...) {
    vn <- deparse(substitute(object)) 
    vn <- gsub(".*\\$", "", vn)
    as.basis(ordered_var(vn, levels = levels(object)))
}


