"shape.statistic" <-
function(tree, norm=NULL) {
	
	if (identical(tree,NULL)) {
		stop("invalid tree","\n")
		
	}	
	if (identical(norm,NULL)) {
		return(sum(log(smaller.clade.spectrum(tree)[,1]-1)))	
	}
	else if (norm=="pda") {
		n<-nrow(tree$merge)
		res<-(sum(log(smaller.clade.spectrum(tree)[,1]-1))-2.03*(n+1)+3.545*sqrt(n))/sqrt(1.570*n*log(n)-5.674*n+3.602*sqrt(n)+14.915)
		return(res)
	}
	else if (norm=="yule") {
		n<-nrow(tree$merge)
		res<-(sum(log(smaller.clade.spectrum(tree)[,1]-1))-1.204*(n+1)+log(n)+2)/sqrt(0.168*(n+1)-0.71)
		return(res)
	}
	else {
		stop("Incorrect argument for 'norm'")
	}
}

