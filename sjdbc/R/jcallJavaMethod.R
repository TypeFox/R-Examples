if (version$language == "R") { 
".JavaMethod" <- 
function(X, METHOD, SIGNATURE, ..., client = FALSE) {
	dotArgs <- list(...)
	# recode Ljava/lang/String to S
	SIGNATURE <- gsub("Ljava/lang/String;", "s", SIGNATURE, fixed=TRUE)
	# recode capital S (short) to T 
	SIGNATURE <- gsub("S", "T", SIGNATURE, fixed=TRUE)
	# recode lowercase s to capital S (String)
	SIGNATURE <- gsub("s", "S", SIGNATURE, fixed=TRUE)
    if (regexpr("^\\((\\[?[STZBCSIJFD])*\\)(\\[?[STZBCSIJFD]|V)$", SIGNATURE)!=1) {
        stop("bad Java method signature")
    }
	
	# tokenize JNI arguments
	sigm <- regexpr("\\(.*\\)", SIGNATURE)
	sigargs <- substring(SIGNATURE, sigm+1, sigm+attr(sigm, "match.length")-2)
	sigm <- gregexpr("\\[?.", sigargs)[[1]]
	sigargs <- substring(sigargs, sigm, sigm+attr(sigm, "match.length")-1)	
	
	# need to wrap array args in .jarray();
	for (i in which(nchar(sigargs)==2)) {
		dotArgs[[i]] <- .jarray(dotArgs[[i]])
	}
	# need to wrap integer args in as.integer explicitly
	for (i in which(sigargs=="I")) {
		dotArgs[[i]] <- as.integer(dotArgs[[i]])
	}
	
	args <- c(list(obj=X,
		   returnSig= gsub(".*)", "", SIGNATURE),
		   method=METHOD), 
		   dotArgs)
	do.call(".jcall", args)
}
}
