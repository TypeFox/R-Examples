part2hi <- function(Part){
    ind = length(Part)
    lastP = Part[[1]]
    if (length(Part[[ind]]) == 0){
        Part[[ind]] = lastP
        for (j in 1:(length(lastP))){
            Part[[ind]][[j]] = rep(1, length(Part[[ind]][[j]]))
        }
    }
    n.levels <- as.vector(sapply(Part[[1]], length))
    n.ll <- n.levels
    dhi = sum(n.levels - 1) + 1
    hi = matrix(0, (ind - 1), dhi)
    nn <- unlist(sapply(1:length(lastP), function(x) rep(names(lastP)[[x]], n.levels[x]-1)))
    colnames(hi) = c("intercept",nn)
    for (i in 2:ind){
        n.levels.new <- sapply(Part[[i]], levels)
        n.levels.new <- as.vector(sapply(n.levels.new, length))
        sum.lev = 0
        for (j in 1:length(n.ll)){
            if (sum(Part[[i]][[j]] != lastP[[j]]) != 0){
                ii = min(which(Part[[i]][[j]] != lastP[[j]]))
                hi[i - 1, sum.lev + ii] = 1
                indeks <- which(Part[[i]][[j]] == Part[[i]][[j]][ii])
                ii2 = min(indeks[indeks != ii])
                if (colnames(hi)[sum.lev + ii2] != colnames(hi)[sum.lev + ii]){
                    hi[i - 1, 1] = -1
                }  else {
                    hi[i - 1, sum.lev + ii2] = -1
                }
            }
            sum.lev = sum.lev + n.ll[j]-1
        }
        n.levels = n.levels.new
        lastP = Part[[i]]
    }
    return(hi)
}
