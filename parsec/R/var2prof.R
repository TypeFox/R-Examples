var2prof <-
function(varmod = lapply(as.list(varlen), function(x) 1:x), varlen = sapply(varmod, length), freq = NULL, labtype = c("profiles", "progressive"), y=NULL) {
        
    profiles <- expand.grid(varmod)
    m <- nrow(profiles)
    
    if(is.null(freq))
        freq <- rep(1, m)
    
    if (length(freq)!=m)
        stop("the length of the vector of frequencies differs from of the number of profiles")
    labtype <- labtype[1]
    if(labtype=="profiles")
        rownames(profiles) <- apply(profiles, 1, function(x) paste(x, collapse=""))
    if(labtype=="progressive")
        rownames(profiles) <- sprintf(paste("P%0", ceiling(log(m, 10)), "i", sep=""), 1:m)
    if(length(y)!=0)
        freq <- apply(profiles, 1, function(z)
            sum(apply(y, 1, function(x) all(z == x)))
        )
    names(freq) <- rownames(profiles)
    res <- list(profiles=profiles, freq=freq)
    class(res) <- "wprof"
    return(res)
}
