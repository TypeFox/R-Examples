intersecZ <- function(Z1, Z2) {
    ord <- rownames(Z1)
    res <- (Z1 * Z2[ord, ord])!=0
    return(res)
}

mrg <- function(lst,
                varmod = lapply(as.list(varlen), function(x) 1:x),
                varlen = sapply(varmod, length)
)
    UseMethod("mrg", object = lst[[1]])

mrg.incidence <- function(lst, varmod = NULL, varlen = NULL)
{  
    n <- length(lst)
    Zres <- lst[[1]]
    if (n > 1) for (i in 2:n)
        Zres <- intersecZ(Zres, lst[[i]])
    class(Zres) <- "incidence"
    
    return(Zres)
}

mrg.character <- function(lst,
    varmod = lapply(as.list(varlen), function(x) 1:x),
    varlen = sapply(varmod, length))
{  
    n <- length(lst)
    Zres <- tmp <- LE2incidence(lst[[1]], varmod)
    if (n > 1) for (i in 2:n) {
        tmp <- LE2incidence(lst[[i]], varmod)
        Zres <- intersecZ(Zres, tmp)
    }
    rm(tmp)
    class(Zres) <- "incidence"
    
    return(Zres)
}