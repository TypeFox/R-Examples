homoteneity <- function (taxa,clustering) 
{
    clustering <- clustify(clustering)
    levels <- levels(clustering)
    clustering <- as.integer(clustering)

    numtyp <- length(table(clustering))
    homo <- rep(NA,numtyp)
    S <- mean(apply(taxa>0,1,sum))
    const <- const(taxa,clustering)
    for (i in 1:numtyp) {
        tmp <- as.numeric(rev(sort(const[,i])))
        homo[i] <- mean(tmp[1:S])
    }
    out <- data.frame(as.character(1:numtyp),homo)
    names(out) <- c('cluster','homoteneity')
    attr(out,'call') <- match.call()
    attr(out,'orig_clustering') <- levels
    out
}
