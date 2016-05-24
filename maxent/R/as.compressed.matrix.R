as.compressed.matrix <- function(DocumentTermMatrix) {
	if (pmatch("matrix.csr",class(DocumentTermMatrix),nomatch=0) > 0) { return(DocumentTermMatrix); }	
	else if (pmatch("TermDocumentMatrix",class(DocumentTermMatrix),nomatch=0) > 0) { DocumentTermMatrix <- t(DocumentTermMatrix); }
    else if (pmatch("DocumentTermMatrix",class(DocumentTermMatrix),nomatch=0) > 0) { flag <- TRUE; }
    else if (pmatch("Matrix",attr(class(DocumentTermMatrix),"package"),nomatch=0) > 0) { return(as.matrix.csr(DocumentTermMatrix)); }
    else if (pmatch("data.frame",class(DocumentTermMatrix),nomatch=0) > 0) { return(as.matrix.csr(as.matrix(DocumentTermMatrix))); }
    else if (pmatch("matrix",class(DocumentTermMatrix),nomatch=0) > 0) { return(as.matrix.csr(DocumentTermMatrix)); }
    else {
    tryCatch(return(as.matrix.csr(as.matrix(DocumentTermMatrix))),error=function(e) stop("Data must be encapsulated using one of the following classes: DocumentTermMatrix or TermDocumentMatrix (package tm), Matrix (package Matrix), matrix.csr (SparseM), data.frame, or matrix"));
    }
    
	ia <- c(1);
	for (n in 1:dim(DocumentTermMatrix)[1]) {
		el <- sum(DocumentTermMatrix$i == n)+ia[length(ia)];
		ia <- append(ia,el);
	}
	
	matrix <- new("matrix.csr",ra=as.numeric(DocumentTermMatrix$v),ja=DocumentTermMatrix$j,ia=as.integer(ia),dimension=dim(DocumentTermMatrix));
	
	return(matrix);
}