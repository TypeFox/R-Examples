classmatch <- function (x,y,type='full')
{
    x <- clustify(x)
    y <- clustify(y)

    tab <- table(x,y)
    total <- sum(tab)
    size <- max(nrow(tab),ncol(tab))

    match <- sum(tab>0)
    pairs <- matrix(0,nrow=match,ncol=3)
    partial <- rep(0,match)
    combo <- rep(0,length(x))
    ord <- matrix(0,nrow=nrow(tab),ncol=ncol(tab))
    running <- 0

    for (i in 1:match) {
        find <- max(tab)
        for (j in 1:nrow(tab)) {
            test <- max(tab[j,])
            if (test == find) {
                col <- which.max(as.vector(tab[j,]))
                pairs[i,] <- c(j,col,tab[j,col])
                tab[j,col] <- 0
                ord[j,col] <- i
                break
            }
        }
    }

    for (i in 1:length(x)) {
        for (j in 1:nrow(pairs)) {
            if (x[i] == pairs[j,1] && y[i] == pairs[j,2]) {
                combo[i] <- j
                break
            }
        }
    }

    partial <- cumsum(pairs[,3])/total
    pairs <- data.frame(pairs)
    names(pairs) <- c('row','column','n')

    if (type=='full')
        res <- list(tab=table(x,y),pairs=pairs,partial=partial,ord=ord,combo=combo)
    else
        res <- list(tab=table(x,y),pairs=pairs,partial=partial[1:size])
    res
}

