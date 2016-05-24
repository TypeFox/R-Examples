# functions for computing percentage tables

# last modified 2014-08-04 by J. Fox

colPercents <- function(tab, digits=1){
    dim <- length(dim(tab))
    if (is.null(dimnames(tab))){
        dims <- dim(tab)
        dimnames(tab) <- lapply(1:dim, function(i) 1:dims[i])
    }
    sums <- apply(tab, 2:dim, sum)
    per <- apply(tab, 1, function(x) x/sums)
    dim(per) <- dim(tab)[c(2:dim,1)]
    per <- aperm(per, c(dim, 1:(dim-1)))
    dimnames(per) <- dimnames(tab)
    per <- round(100*per, digits)
    result <- abind(per, Total=apply(per, 2:dim, sum), Count=sums, along=1)
    names(dimnames(result)) <- names(dimnames(tab))
    result
}

rowPercents <- function(tab, digits=1){
    dim <- length(dim(tab))
    if (dim == 2) return(t(colPercents(t(tab), digits=digits)))
    tab <- aperm(tab, c(2,1,3:dim))
    aperm(colPercents(tab, digits=digits), c(2,1,3:dim))
}

totPercents <- function(tab, digits=1){
    dim <- length(dim(tab))
    if (is.null(dimnames(tab))){
        dims <- dim(tab)
        dimnames(tab) <- lapply(1:dim, function(i) 1:dims[i])
    }
    tab <- 100*tab/sum(tab)
    tab <- cbind(tab, rowSums(tab))
    tab <- rbind(tab, colSums(tab))
    rownames(tab)[nrow(tab)] <- "Total"
    colnames(tab)[ncol(tab)] <- "Total"
    round(tab, digits=digits)
}