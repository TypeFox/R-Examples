monophyly <- 
function (phy, sppVector, pp = NA, singletonsMono = TRUE) 
{
    res <- list()
    xxx <- lapply(unique(sppVector), function(y) which(sppVector == 
        y))
    sppTab <- sapply(xxx, length)
    singletons <- which(sppTab == 1)
    nonSingletons <- which(sppTab != 1)
    ifelse(is.na(pp), yyy <- prop.part(phy), yyy <- pp)
    zzz <- sapply(yyy, length)
    defNon <- which(!sppTab %in% zzz)
    poss <- which(sppTab %in% zzz)
    for (i in poss) {
        res[i] <- NA
        for (j in 1:length(yyy[which(zzz == sppTab[i])])) res[[i]][j] <- sum(as.numeric(!xxx[[i]] %in% 
            yyy[which(zzz == sppTab[i])][[j]]))
    }
    out <- sapply(res, function(x) as.logical(sum(as.numeric(x < 1))))
    if(is.list(out)) out <- rep(singletonsMono, length(singletons))
    out[defNon] <- FALSE
    out[singletons] <- singletonsMono
    out
}


