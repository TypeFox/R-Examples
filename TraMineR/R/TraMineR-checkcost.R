
TraMineR.checkcost <- function(sma, seqdata, with.missing, indel, tol = NULL) {
	if(is.null(tol)){
		if(!missing(indel)){
			tol <- 1e-7*indel
		}else{
			tol <- 1e-7
		}
	}
	alphabet <- attr(seqdata,"alphabet")
	## Gaps in sequences
	if (with.missing) {
		alphabet <- c(alphabet,attr(seqdata,"nr"))
	}
	alphsize <- length(alphabet)
	
	if(length(dim(sma))==2){
		dim(sma) <- c(dim(sma),1)
	}
	else {
		if(ncol(seqdata)!=dim(sma)[3]){
			stop(" [!] size of substitution cost matrix must be ", alphsize,"x", alphsize, "x", ncol(seqdata), call. = FALSE)
		}
	}
	if (!missing(indel) && indel <= 0) {
		stop(" [!] indel cost should be positive")
	}
	for(i in 1:dim(sma)[3]) {
		sm <- sma[,,i]
		if (nrow(sm)!=alphsize | ncol(sm)!=alphsize) {
			stop(" [!] size of substitution cost matrix must be ", alphsize,"x", alphsize, call. = FALSE)
		}
		if (any(sm<0)) {
			stop(" [!] Negative substitution costs are not allowed", call. = FALSE)
		}
		
		if (any(diag(sm)!=0)) {
			stop(" [!] All element on the diagonal of sm (substitution cost) should be equal to zero")
		}
		triangleineq <- checktriangleineq(sm, warn=FALSE, indices=TRUE, tol=tol)
		## triangleineq contain a vector of problematic indices.
		if (!is.logical(triangleineq)) {
			
			warning(" [!] at least, one substitution cost doesn't respect the triangle inequality.\n",
        			" [!] replacing ", alphabet[triangleineq[1]], " with ", alphabet[triangleineq[3]], " (cost=", format(sm[triangleineq[1], triangleineq[3]]), 
					") and then ", alphabet[triangleineq[3]], " with ", alphabet[triangleineq[2]], " (cost=", format(sm[triangleineq[3], triangleineq[2]]), 
					")\n [!] costs less than replacing directly ", alphabet[triangleineq[1]], " with ", alphabet[triangleineq[2]], 
					" (cost=", format(sm[triangleineq[1], triangleineq[2]]), ")\n",
					" [!] total difference ([", alphabet[triangleineq[1]], "=>", alphabet[triangleineq[3]],
					"] + [", alphabet[triangleineq[3]], "=>", alphabet[triangleineq[2]],"] - [", alphabet[triangleineq[1]], "=>", alphabet[triangleineq[2]],"]): ",
					format(sm[triangleineq[1], triangleineq[3]]+sm[triangleineq[3], triangleineq[2]]-sm[triangleineq[1], triangleineq[2]]), call. = FALSE)
		}
		if(!missing(indel) && indel<=0){
			stop(" [!] indel should be greater than zero")
		}
		if (!missing(indel) && any(sm>2*indel)) {
			warning("Some substitution cost are greater that two times the indel cost.",
				" Such substitution cost will thus never be used.", call. = FALSE)
		}
		
		## Testing for symmetric matrix
		if (any((sm-t(sm))>tol)) {
			warning("The substitution cost matrix is not symmetric.", call. = FALSE)
		}
	}
}