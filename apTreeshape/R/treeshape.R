"treeshape" <-
function(nodes, names){

#nodes est une matrice a deux colonnes
	if (!(is.matrix(nodes)) || (ncol(nodes) != 2)) {
		stop("invalid argument")
	}
#names est un vecteur a n elements.
	if (missing(names)) {
		names=paste("tip",as.character(1:(nrow(nodes)+1)))
	}
	
	if (length(names)!=nrow(nodes)+1){
		stop("Wrong number of elements for vector names")
	}
	
	res <- list(merge=nodes, names=names)
	class(res)<-'treeshape'
	
	res
}

