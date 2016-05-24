is.subseqelist<-function(s){
#   return(.Call("istmrsequence",s, PACKAGE="TraMineR"))
	return(inherits(s, "subseqelist"))
}

createsubseqelist<-function(seqe, constraint, subseq, data, type="frequent"){
	if(!is.seqelist(seqe)) {
		stop(" [!] seqe should be a seqelist")
	}
	ret <- list()
	ret$seqe <- seqe
	ret$constraint <- constraint
	ret$subseq <- subseq
	ret$data <- as.data.frame(data, optional=TRUE)
	ret$type <- type
	class(ret$subseq) <- c("seqelist", "list")
#  attr(ret$subseq,"dictionnary")<-attr(seq,"dictionnary")
	class(ret) <- "subseqelist"
	return(ret)
}
print.subseqelist<-function(x,...){
	z <- data.frame(data.frame(Subsequence=as.character(x$subseq), check.names=FALSE), x$data, row.names=NULL, check.names=FALSE)
	print(z, ...)
	cat("\nComputed on", length(x$seqe), "event sequences\n")
	print(x$constraint,...)
}

"[.subseqelist" <- function(x, i,j, drop=FALSE) {
    # If only 1 subscript is given, the result will still be a Surv object
    #  If the second is given extract the relevant columns as a matrix
	if (missing(j)) {
		ret <- createsubseqelist(x$seqe, x$constraint, x$subseq[i,drop=drop], data=x$data[i,], type=x$type)
		if(!is.null(x$labels)) {
			ret$labels<-x$labels
		}
		if(!is.null(x$bonferroni)) {
			ret$bonferroni<-x$bonferroni
		}
		class(ret)<-class(x)
		return(ret)
	} else {
		class(x) <- NULL
		NextMethod("[")
	}
}

plot.subseqelist<-function(x, freq=NULL, cex=1,...){
	slegend <- as.character(x$subseq)
	if(is.null(freq)) {
		freq<-x$data[,1]
	}
	barpos <- barplot(freq, names.arg=c(""), ...)
	text(x=barpos, y=0.02, labels=slegend, srt=90, adj=c(0,0.5), cex=cex)
}