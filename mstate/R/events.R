`events` <- function(msdata)
{
    trans <- attr(msdata, "trans")
    K <- nrow(trans)
    absorbing <- which(apply(is.na(trans), 1, all))
    if (!is.null(dimnames(trans))) states <- dimnames(trans)[[1]]
    else states <- as.character(1:K)
    from <- factor(msdata$from,levels=1:K,labels=states)
    to <- factor(msdata$to,levels=1:K,labels=states)
    tbl <- table(from[msdata$status==1],to[msdata$status==1],dnn=c("from","to"))
    counts <- tbl
    tbl <- table(from,to,dnn=c("from","to"))
    total <- apply(tbl,1,max)
    total[absorbing] <- apply(counts[,absorbing,drop=FALSE], 2, sum)
    noevent <- total - apply(counts,1,sum)
    counts <- cbind(counts,noevent,total)
    dn <- dimnames(counts)
    dn[[2]][(K+1):(K+2)] <- c("no event","total entering")
    names(dn) <- c("from","to")
    dimnames(counts) <- dn
    class(counts) <- "table"
    freqs <- (counts/total)[,-(K+2)]
    class(freqs) <- "table"
    return(list(Frequencies=counts,Proportions=freqs))
}
