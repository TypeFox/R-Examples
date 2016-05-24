cem.match <- function (data, verbose = 0) 
{
    vnames <- colnames(data)
    n <- dim(data)[1]
    nv <- dim(data)[2]
    if (verbose > 1) {
        cat("\nmatching on variables:")
        cat(paste(vnames))
        cat("\n")
    }
    xx <- apply(data, 1, function(x) paste(x, collapse = "\r"))
    tab <- table(xx)
    st <- names(tab)
    strata <- match(xx,st)
    n.strata <- length(st)
    return(invisible(list(call = match.call(), strata = strata, 
        n.strata = n.strata, vars = vnames)))   
}
