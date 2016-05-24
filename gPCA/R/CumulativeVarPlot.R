CumulativeVarPlot <-
function(out,ug="unguided",...){
	# out - output of gPCA.batchdetect
	# ug - "guided" or "unguided" principal components plots?

	if (ug=="unguided") {
		plot(out$cumulative.var.u,type="l",xlab="Principal Component",ylab="Cumulative Variance",ylim=c(0,1),...)
	} else {
		plot(out$cumulative.var.g,type="l",xlab="Principal Component",ylab="Cumulative Variance",ylim=c(0,1),xaxt="n",...)
		axis(1,1:out$b,paste(1:out$b))
	}
	
}
