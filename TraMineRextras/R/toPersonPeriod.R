toPersonPeriod <- function(seqdata){
    names(seqdata) <- paste(names(seqdata),"..", sep="")
	ids <- 1:nrow(seqdata)
	if(!is.null(rownames(seqdata))){
		ids <- rownames(seqdata)
	}
    pp <- data.frame(id=rep(ids, ncol(seqdata)), state=unlist(seqdata), 
					timestamp =sort(rep(0:(ncol(seqdata)-1), nrow(seqdata))))
	if(!is.null(attr(seqdata, "void"))){
		pp <- pp[pp$state != attr(seqdata, "void"), ]
	}
	if(!is.null(attr(seqdata, "nr"))){
		pp$state[pp$state==attr(seqdata, "nr")] <- NA
	}
    pp <- pp[order(pp$id, pp$timestamp, pp$state), ]
	return(pp)
}
