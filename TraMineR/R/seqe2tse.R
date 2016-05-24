seqe2TSE <- function(seqe){
	tse <- .Call("tmrseqetotse",  unlist(list(seqe)), PACKAGE="TraMineR")
	ll <- levels(seqe)
	tse <- data.frame(id=tse[[1]], timestamp=tse[[2]], event=factor(tse[[3]], levels=1:length(ll), labels=ll))
	return(tse)
}