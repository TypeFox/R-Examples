
summary.stslist <- function(object,...) {

	alphabet <- alphabet(object)
	nbstates <- length(alphabet)
	cpal <- cpal(object)
	labels <- attr(object,"labels")
	nr <- attr(object,"nr")
	void <- attr(object,"void")
	weights <- attr(object, "weights")
	TraMineR.version <- attr(object, "Version")

	nbseq <- seqdim(object)[1]
	seql <- seqlength(object)
	nuseq <- nrow(unique(object))

	if (!is.null(TraMineR.version)) {
		cat(" [>] sequence object created with TraMineR version",TraMineR.version,"\n")
	}
	cat(" [>]", nbseq, "sequences in the data set,",nuseq, "unique","\n")

	## weights
	if (!is.null(weights) && !all(weights==1)) {
		cat(" [>] sum of weights: ", round(sum(weights),2), " - min/max: ",
			min(weights),"/",max(weights),"\n", sep="")
	}
	cat(" [>] min/max sequence length: ",min(seql),"/",max(seql),"\n", sep="")

	## Alphabet
	cat(" [>] alphabet (state labels): ","\n")
	maxstatedisplay <- 12
	for (i in 1:min(nbstates,maxstatedisplay))
		cat("     ",i, "=", alphabet[i], " (", labels[i], ")","\n", sep="")
	if (nbstates>12) message("      ...")
	cat(" [>] dimensionality of the sequence space:", (nbstates-1)*max(seql),"\n")
	cat(" [>] colors:", paste(1:nbstates,cpal,collapse=" ",sep="="),"\n")

	if (any(object==nr)) {	
		cat(" [>] symbol for missing state:",nr,"\n")
	}

	if (any(object==void)) {	
		cat(" [>] symbol for void element:",void,"\n")
	}
}

