seqelatex <- function(seqe){
	xx <- gsub("-([0-9.]+)-?", " \\\\xrightarrow{\\1} ", as.character(seqe))
	xx <- gsub("-", " \\\\rightarrow ", xx)
	xx <- gsub("\\(([^\\)]+)\\)", "(\\\\mbox{\\1})", xx)
	return(paste("$",xx, "$", sep=""))
}

seqelatex.small <- function(seqe){
	xx <- gsub("-([0-9.]+)-?", " \\\\xrightarrow{\\1} ", as.character(seqe))
	xx <- gsub("-", " \\\\rightarrow ", xx)
	xx <- gsub("\\(([^\\)]+)\\)", "(\\\\mbox{\\1})", xx)
	return(paste("{\\small $",xx, "$}", sep=""))
}

seqelatex.script <- function(seqe){
	xx <- gsub("-([0-9.]+)-?", " \\\\xrightarrow{\\1} ", as.character(seqe))
	xx <- gsub("-", " \\\\rightarrow ", xx)
	xx <- gsub("\\(([^\\)]+)\\)", "(\\\\mbox{\\1})", xx)
	return(paste("{\\scriptsize $",xx, "$}", sep=""))
}


##plotsubseqelist<-function(x, freq=NULL, cex=1,file="tplot.tex", standAlone =TRUE, ...){
##	tikz(file=file,  standAlone =  standAlone)
##	slegend <- seqelatex(x$subseq)
##	if(is.null(freq)) {
##		freq<-x$data[,1]
##	}
##	barpos <- barplot(freq, names.arg=c(""), ...)
##	text(x=barpos, y=0.02, labels=slegend, srt=90, adj=c(0,0.5), cex=cex)
##	dev.off()
##}
