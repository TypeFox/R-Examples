prolonged.sed.bouts.min <-
function(posture,epoch=1,n) {	
	acts <- posture == 0
	lengths.of.continuous.bouts.of.sed <- apply(as.matrix(strsplit(paste(acts, collapse=""), split="FALSE", fixed=TRUE)[[1]]), 1, function(x) {nchar(x)/4})
	lengths <- lengths.of.continuous.bouts.of.sed > n*(60/epoch)
	return(sum(lengths.of.continuous.bouts.of.sed[lengths])/(60/epoch))
	}

