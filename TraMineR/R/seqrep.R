## =============================================
## Representative sequence of a set of sequences
## =============================================

seqrep <- function(seqdata, criterion="density", score=NULL, decreasing=TRUE, 
	trep=0.25, nrep=NULL, tsim=0.10, 
	dmax=NULL, dist.matrix=NULL, weighted=TRUE, ...) {

	if (!inherits(seqdata,"stslist")) 
		stop("data is not a sequence object, see seqdef function to create one", call.=FALSE)

	slength <- ncol(seqdata)
	statelist <- alphabet(seqdata)

	weights <- attr(seqdata, "weights")

	if (!weighted || is.null(weights)) { weights <- rep(1, nrow(seqdata)) }
	if (all(weights==1)) { weighted <- FALSE }

	## Distance matrix
	if (missing(dist.matrix) || is.null(dist.matrix))
		dist.matrix <- seqdist(seqdata, ...)


	if (is.null(score)) {
		## ====================
		## Max representativity
		## ====================
		if (criterion=="mscore") {
			message(" [!] criterion still under development")

			## State distribution
			freq <- seqstatd(seqdata)$Frequencies

			score <- apply(seqdata,1, TraMineR.mscore, slength, statelist, freq)
			decreasing <- TRUE
		} 
		## ===============
		## Max probability
		## ===============
		else if (criterion=="prob") {

			score <- seqlogp(seqdata)
			decreasing <- FALSE
		}
	}

	## Getting the representatives
	rep <- dissrep(dist.matrix, criterion=criterion, score=score, 
		decreasing=decreasing, trep=trep, nrep=nrep, tsim=tsim, dmax=dmax, weights=weights)
	
	## Occurence of the representative sequence
	nds <- nrow(unique(seqdata))
	message(" [>] ", nds, " distinct sequence(s)")

	## ============
	## Final object
	## ============
	res <- seqdata[rep,]
	rownames(res) <- paste("[",1:nrow(res),"]", sep="")
	class(res) <- c("stslist.rep", class(res))

	attr(res, "nbseq") <- attr(rep, "n")
	attr(res, "criterion") <- criterion
	attr(res, "dmax") <- attr(rep,"dmax")
	attr(res, "Index") <- as.vector(rep)
	attr(res, "Scores") <- attr(rep,"Scores")
	attr(res, "Distances") <- attr(rep,"Distances")
	attr(res, "Statistics") <- attr(rep,"Statistics")
	attr(res, "Quality") <- attr(rep,"Quality") 
	attr(res, "weighted") <- weighted

	return(res)
}	
	
