getlambda <- function(...) {
    rels <- as.list(match.call())[-1]
    rels <- lapply(rels, as.character)
    rels <- lapply(rels, function(x) {
        if(length(x) == 1)
            x <- c("=", x, x)
        return(x)
    })
    rels <- matrix(unlist(rels), nrow=3)
    vars <- unique(as.vector(rels[-1,]))
    n <- length(vars)
    res <- matrix(FALSE, n, n, dimnames = list(vars, vars))
    apply(rels, 2, function(x){
        if (x[1] == "<")
            res[x[2], x[3]] <<- TRUE
        if (x[1] == ">")
            res[x[3], x[2]] <<- TRUE
    })
    class(res) <- "cover"
    res <- cover2incidence(res)
    if (!is.partialorder(res))
        stop("the introduced relations do not generate a partial order")
return(res)
}