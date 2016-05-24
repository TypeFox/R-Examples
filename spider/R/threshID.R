threshID <- function(distobj, sppVector, threshold = 0.01){
	distobj <- as.matrix(distobj)
	diag(distobj) <- NA
	output <- rep(NA, length(sppVector))
	aa <- apply(distobj, MARGIN=2, FUN=function(x) which(x < threshold))
	if(is.matrix(aa)) aa <- c(unname(as.data.frame(aa))) else aa
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

######################
#Return the name of the result

#~ idBOLD <- function(distobj, sppVector, threshold = 0.01, highThreshold = 0.05){
	#~ distobj <- as.matrix(distobj)
	#~ diag(distobj) <- NA
	#~ output <- rep(NA, length(sppVector))
	#~ for(i in 1:length(sppVector)){
		#~ high <- sppVector[which(distobj[,i] < highThreshold)]
		#~ low <- sppVector[which(distobj[,i] < threshold)]

		#~ if(length(high) == 0) output[i] <- NA else {
		  #~ if(length(unique(low)) == 1) output[i] <- unique(low) else output[i] <- FALSE}
	#~ }
	#~ output
#~ }
