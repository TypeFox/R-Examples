## compute pairwise fold changes
pairwise.fc <- function (x, g, ave = mean, log = TRUE, base = 2, 
                         mod.fc = TRUE, ...) {
    g <- factor(g)
    compare.levels <- function(i, j) {
        xi <- x[as.integer(g) == i]
        xj <- x[as.integer(g) == j]
        if(log){
            logFC <- ave(xj, ...) - ave(xi, ...)
            FC <- base^logFC
        }else{
            FC <- ave(xj, ...)/ave(xi, ...)
        }
        if(mod.fc)
          return(ifelse(FC > 1, FC, -1/FC))
        else
          return(FC)
    }
    ix <- seq_along(levels(g))
    names(ix) <- levels(g)
    pp <- outer(ix[-1], ix[-length(ix)], function(ivec, jvec) sapply(seq_along(ivec), 
        function(k) {
            i <- ivec[k]
            j <- jvec[k]
            if (i > j) 
                compare.levels(i, j)
            else NA
        }))
    pp <- pp[lower.tri(pp, diag = TRUE)]

    nr <- length(levels(g))
    Names <- numeric(choose(nr, 2))
    Levels <- levels(g)
    count <- 0
    for(i in seq_len(nr-1))
    for(j in (i+1):nr){
        count <- count + 1
        Names[count] <- paste(Levels[i], Levels[j], sep = " vs ")
    }

    names(pp) <- Names
    return(pp)
}
