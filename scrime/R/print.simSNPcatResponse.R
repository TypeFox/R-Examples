print.simSNPcatResponse<-function(x,justify=c("left","right"),spaces=2,...){
	if(spaces<0)
		stop("spaces must be non-negative.")
	tab<-x$tab.explain
	justify<-match.arg(justify)
	tab<-format(tab,justify=justify)
	if(justify=="left"){
		n<-nchar(tab$IA[1])-2
		if(spaces>n)
			stop("spaces has to be at most ",n,".")
		name.ia<-paste(paste(rep(" ",spaces),collapse=""),"IA",
			paste(rep(" ",n-spaces),collapse=""),sep="")
		colnames(tab)[ncol(tab)]<-name.ia
	}
	print(tab)
}

