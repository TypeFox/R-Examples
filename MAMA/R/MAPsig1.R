MAPsig1<-function(unique.pat,value.dis, iter=1000) 
{
cat("Permutation: \n")
    n <- apply(value.dis, 2, sum)
    n.soft <- patternMatch(value.dis, unique.pat)
    n.strong <- patternMatch.strong(value.dis, unique.pat)
    n.pat <- length(unique.pat)
    res.random <- matrix(0, length(unique.pat), iter)
    res.random.strong <- matrix(0, length(unique.pat), iter)
    n.entity <- dim(value.dis)[2]
    genes <- 1:nrow(value.dis)
    for (l in 1:iter) {
        if (l %% 50 == 0) cat(l, "\n")
        Sgenex.random <- value.dis * 0
        for (i in 1:n.entity) {
            Sgenex.random[sample(genes, n[i]), i] <- 1
        }
        res.random[, l] <- patternMatch(Sgenex.random, unique.pat)
        res.random.strong[, l] <- patternMatch.strong(Sgenex.random, 
            unique.pat)
    }
    rownames(res.random) <- unique.pat
    rownames(res.random.strong) <- unique.pat
    p.soft <- array(0, n.pat)
    p.strong <- array(0, n.pat)
    for (i in 1:n.pat) {
        p.soft[i] <- length(which(res.random[i, ] >= n.soft[i]))/iter
        p.strong[i] <- length(which(res.random.strong[i, ] >= 
            n.strong[i]))/iter
    }
    res <- data.frame(p.soft, p.strong)
    return(res)

}
