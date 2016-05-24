getMatFuzzy <-
function(score, probs, check=TRUE){
	n.mat <- length(probs)
	if(!length(probs) %in% 2:3)
		stop("probs must consist of either 2 or 3 elements.")
	if(is.character(score)){
		if(length(score) != 1)
			stop("If score is of class character, it must contain only of one element.")
		type <- match.arg(score, c("additive", "dominant", "recessive"))
		score <- switch(type, additive=0:2, dominant=c(0,1,1), recessive=c(0,0,1))
		if(n.mat==2)
			score <- score[-1]
	}
	else{
		if(!length(score) %in% 2:3)
			stop("score must consist either of 2 or 3 elements.")
		if(!is.numeric(score))
			stop("score must be a numeric vector.")
		if(any(score<0))
			stop("All scores must be non-negative.")
		if(length(score) == 3 && score[1] != 0)
			stop("If score contains three values, the first one must be zero.")
	}
	if(n.mat != length(score))
		stop("score and probs must consist of the same number of elements.")
	if(!is.list(probs))
		stop("probs must be a list.")
	if(check){
		tmp <- 0
		d1 <- dim(probs[[1]])
		rn <- rownames(probs[[1]])
		cn <- colnames(probs[[1]])
		for(i in 1:n.mat){
			if(!is.matrix(probs[[i]]))
				stop("All objects in probs must be matrices.")
			if(any(dim(probs[[i]]) != d1))
				stop("The dimensions of the matrices in probs are not identical.")
			if(!is.null(rn) && any(rownames(probs[[i]]) != rn))
				stop("The row names of the matrices in probs differ between the matrices.")
			if(!is.null(cn) && any(colnames(probs[[i]]) != cn))
				stop("The column names of the matrices in probs differ between the matrices.")
			if(any(probs[[i]]<0) || any(probs[[i]]>1))
				stop("The values in the matrices in probs must be between 0 and 1.")
			tmp <- tmp + probs[[i]]
		}
	}
	if(n.mat==2){
		if(check && any(tmp>1))
			stop("Some of the sums over the probabilities in the matrices in probs are larger than 1.")
		return(score[1] * probs[[1]] + score[2] * probs[[2]])
	}
	if(check && any(round(tmp,6)!=1))
		stop("The probabilities in the matrices must always sum up to 1.")
	score[2] * probs[[2]] + score[3] * probs[[3]]
}

