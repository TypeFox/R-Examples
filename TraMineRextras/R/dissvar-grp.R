dissvar.grp <- function(mdis, group=NULL, ...){

	isdist <- inherits(mdis, "dist")
	if (isdist) {
		n <- attr(mdis, "Size")
	} else if (is.matrix(mdis)) {
		n <- nrow(mdis)
	} else {
		stop("mdis argument should be a dist object or a dissimilarity matrix")
	}

    grp <- group
    if (is.null(grp)) {
        grp <- rep(1, n)
    }
    
    if (length(grp) != n){
        stop("length(grp) not compatible with size of mdis",
            call. = FALSE)
        }

    levg <- levels(grp <- factor(grp))

    v <- vector("double",length(levg))
    ## We need to transform into a matrix to subset
    ## an alternative would be use the subset.dist function from package cba
    if (isdist) mdis <- as.matrix(mdis)
    for (i in 1:length(levg))
        {
        ig <- which(grp==levg[i])
        v[[i]] <- dissvar(mdis[ig,ig], ...)
        }
    names(v) <- levg
    return(v)
}
