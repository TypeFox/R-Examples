findwtdinteraction.mira <- function(reg, across, by=NULL, at=NULL, acrosslevs=NULL, bylevs=NULL, atlevs=NULL, weight=NULL, dvname=NULL, acclevnames=NULL, bylevnames=NULL, atlevnames=NULL, stdzacross=FALSE, stdzby=FALSE, stdzat=FALSE, limitlevs=20, type="response"){
    predset <- lapply(reg$analyses, function(g) findwtdinteraction.default(g, across, by, at, acrosslevs=acrosslevs, bylevs=bylevs, atlevs=atlevs, weight=weight, dvname=dvname, bylevnames=bylevnames, atlevnames=atlevnames, acclevnames=acclevnames, stdzacross, stdzby, stdzat, limitlevs=limitlevs))
    allmns <- lapply(predset, function(m) m$Means)
    if(table(table(sapply(allmns, length)))!=1)
        stop("at variable values are inconsistent across imputations, please set atlevs before running")
    if(table(table(unlist(sapply(allmns, function(x) sapply(x, function(y) dim(y)[1])))))!=1)
        stop("by variable values are inconsistent across imputations, please set bylevs before running")
    if(table(table(unlist(sapply(allmns, function(x) sapply(x, function(y) dim(y)[2])))))!=1)
        stop("across variable values are inconsistent across imputations, please set acrosslevs before running")
    allses <- lapply(predset, function(m) m$SEs)
    allresp <- sapply(predset, function(m) m$Resp)
    imputations <- length(allmns)
    nlat <- length(allmns[[1]])
    nlby <- dim(allmns[[1]][[1]])[1]
    nlacross <- dim(allmns[[1]][[1]])[2]
    impmns <- lapply(1:nlat, function(a) sapply(1:nlacross, function(c) sapply(1:nlby, function(b) mean(sapply(1:imputations, function(i) allmns[[i]][[a]][b,c])))))
    impses <- lapply(1:nlat, function(a) sapply(1:nlacross, function(c) as.numeric(sapply(1:nlby, function(b) pool.scalar(sapply(1:imputations, function(i) allmns[[i]][[a]][b,c]), sapply(1:imputations, function(i) allses[[i]][[a]][b,c]))["t"]))))
    for(i in 1:length(impmns)){
        if(!is.vector(impmns[[i]])){
            colnames(impmns[[i]]) <- colnames(impses[[i]]) <- colnames(allmns[[1]][[1]])
            rownames(impmns[[i]]) <- rownames(impses[[i]]) <- rownames(allmns[[1]][[1]])
        }
        else
            names(impmns[[i]]) <- names(impses[[i]]) <- colnames(allmns[[1]][[1]])
    }
    names(impses) <- names(impmns) <- names(allmns[[1]])
    out <- NULL
    out$RespMns <- sapply(as.data.frame(t(allresp)), function(x) try(mean(as.numeric(x), na.rm=TRUE)))
    out$Resp <- allresp
    out$Meta <- predset[[1]]$Meta
    out$Means <- impmns
    out$SEs <- impses
    class(out) <- "interactpreds"
    out
}
