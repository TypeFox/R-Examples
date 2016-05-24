monophylyBoot <- 
function (phy, sppVector, DNAbin, thresh = 0.7, reroot = TRUE, pp = NA, singletonsMono = TRUE, reps = 1000, block = 3) 
{
    res <- list()
    xxx <- lapply(unique(sppVector), function(y) which(sppVector == 
        y))
    if(reroot){
	testTr <- nj(dist.dna(DNAbin))
	maxInt <- max(testTr$edge.length[testTr$edge[,2] > length(testTr$tip.label)])
	nodeRoot <- testTr$edge[which(testTr$edge.length == maxInt), 2]
	boot <- boot.phylo(phy, DNAbin, function(x) root(nj(dist.dna(x)), node = nodeRoot, resolve.root=TRUE), B = reps, block = block)/reps
	} else boot <- boot.phylo(phy, DNAbin, function(x) nj(dist.dna(x)), B = reps, block = block)/reps
    sppTab <- sapply(xxx, length)
    singletons <- which(sppTab == 1)
    nonSingletons <- which(sppTab != 1)
    ifelse(is.na(pp), yyy <- prop.part(phy), yyy <- pp)
    zzz <- sapply(yyy, length)
    defNon <- which(!sppTab %in% zzz)
    poss <- which(sppTab %in% zzz)
    bb <- lapply(sppTab, function(x) boot[which(zzz == x)])
    tt <- lapply(sppTab, function(x) which(zzz == x))
    for(i in poss){
	res[i] <- NA
	for(j in 1:length(tt[[i]])){
		res[[i]][j] <- sum(as.numeric(!xxx[[i]] %in% yyy[[ tt[[i]][j] ]]))
		}
	}
    bc <- sapply(res, function(x) which(x == 0))
    bootCheck <- sapply(1:length(bc), function(x) bb[[x]][bc[[x]][1]])
    out <- sapply(res, function(x) as.logical(sum(as.numeric(x < 
        1))))
    if(is.list(out)) out <- rep(singletonsMono, length(singletons))
    out[defNon] <- FALSE
    out[singletons] <- singletonsMono
    out[bootCheck < thresh] <- FALSE
    store <- list(results= out, BSvalues = boot)
    print(out)
    invisible(store)
}
