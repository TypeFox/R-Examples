"crnames" <- 
function(dnames,idx) {
	rownams <- dnames[[1]]
	colnams <- dnames[[2]]
	if (is.logical(idx)) idx <- which(idx)
	nrw <- length(rownams)
	rnm <- idx %% nrw
	rnm <- replace(rnm,(rnm==0),nrw)
	cnm <- ceiling(idx/nrw)
	out <- matrix(rep(NA,2*length(idx)),ncol=2)
	out[,1] <- rownams[rnm]
	out[,2] <- colnams[cnm]
	colnames(out) <- c("row.name","col.name")
	out
}

