print.farules <- function(x, ...) {
    if (!is.farules(x)) {
        stop("'x' must be an instance of class 'farules'")
    }
    if (length(x$rules) <= 0) {
        cat('Empty rule base\n')
    } else {
        r <- x$statistics
        rownames(r) <- sapply(x$rules, function(rule) {
            ante <- rule[-1]
            conseq <- rule[1]
            paste(paste(ante, collapse=' & '),
                conseq, 
                sep=' => ')
        })
        print(r)
    }
}

