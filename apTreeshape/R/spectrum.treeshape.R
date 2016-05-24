"spectrum.treeshape" <-
function(tree) {
	if (identical(tree,NULL)) {
		stop("invalid tree","\n")
	}

	tmp=smaller.clade.spectrum(tree)
	tmp=tmp[,1]
	res=c()
	for (i in 2:(nrow(tree$merge)+1)){
		res=c(sum(tmp==i), res)
	}
	res
}

