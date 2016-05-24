bestCloseMatch <- function(distobj, sppVector, threshold = 0.01){
	distobj <- as.matrix(distobj)
	diag(distobj) <- NA
	output <- rep(NA, length(sppVector))
	aa <- apply(distobj, MARGIN=2, FUN=function(x) which(x == min(x, na.rm = TRUE)))
	bb <- lapply(aa, function(x) unique(sppVector[x]))
	cc <- sppVector == bb
	dd <- sapply(1:length(sppVector), function(x) sppVector[x] %in% bb[[x]])
	ee <- apply(distobj, MARGIN=2, FUN=function(x) min(x, na.rm = TRUE))
	output[which(cc & dd)] <- "correct"
	output[which(!cc & !dd)] <- "incorrect"
	output[which(!cc & dd)] <- "ambiguous"
	output[which(ee > threshold)] <- "no id"
	output
}