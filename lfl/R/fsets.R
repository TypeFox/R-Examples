fsets <- function(x, vars, specs) {
    if (!is.matrix(x) || !is.numeric(x)) {
        stop("'x' must be a numeric matrix")
    }
    if (!is.vector(vars) || ncol(x) != length(vars) || any(colnames(x) != names(vars))) {
        stop("'vars' must be a vector with names equal to 'colnames(x)'")
    }
    if (!is.matrix(specs) || !is.numeric(specs) || ncol(x) != ncol(specs) 
            || ncol(x) != nrow(specs) || any(colnames(x) != colnames(specs))
            || any(colnames(x) != rownames(specs))) {
        stop("'specs' must be a numeric matrix with colnames and rownames equal to 'colnames(x)'")
    }
    return(structure(x, 
                     class=c('fsets', class(x)),
                     vars=vars,
                     specs=specs))
}
