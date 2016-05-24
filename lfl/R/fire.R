fire <- function(x, 
                 rules,
                 tnorm=c("goedel", "goguen", "lukasiewicz"),
                 onlyAnte=TRUE,
                 parallel=FALSE) {
    if (is.vector(x)) {
        x <- matrix(x, nrow=1, dimnames=list(NULL, names(x)))
    }
    if (!is.matrix(x) || !is.numeric(x)) {
        stop("'x' must be a numeric vector or matrix")
    }
    if (min(x) < 0 || max(x) > 1) {
        stop("Values of 'x' must be truth values in the interval [0,1]")
    }
    if (is.farules(rules)) {
        rules <- rules$rules
    } else if (is.vector(rules) && is.character(rules)) {
        rules <- list(rules)
    }
    if (!is.list(rules)) {
        stop("'rules' must be a list of rules")
    }
    if (!is.function(tnorm)) {
        tnorm <- match.arg(tnorm)
        tnorm <- .tnorms[[tnorm]]
    }

    if (onlyAnte) {
        rules <- antecedents(rules)
    }

    innerLoop <- function(rule) {
        apply(subset(x, select=rule), 1, tnorm)
    }

    rule <- NULL
    if (parallel) {
        res <- foreach(rule=rules) %dopar% { innerLoop(rule) }
    } else {
        res <- foreach(rule=rules) %do% { innerLoop(rule) }
    }

    return(res)
}
