subset.linkdat <- function(x, subset=x$orig.ids, ..., markers=seq_len(x$nMark)) {
	x = removeMarkers(x, setdiff(seq_len(x$nMark), markers))
	xframe = as.matrix(x)
	
	newfr = xframe[xframe[, 'ID'] %in% subset, ,drop=F]
	newfr[!(newfr[, 'FID'] %in% subset), 'FID'] = 0  # set FID=0 if father is not in subset
	newfr[!(newfr[, 'MID'] %in% subset), 'MID'] = 0  # set MID=0 if mother is not in subset
 
	restore_linkdat(newfr, attributes(xframe))
}
