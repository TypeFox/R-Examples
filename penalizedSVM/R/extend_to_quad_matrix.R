
.extend.to.quad.matrix <-
function (mat, labels.universe=NULL){
	# extend to quadratic matrix
	# mat - matrix
	# labels.universe - all labels to be appeared in the quadratic matrix
	
	# if label's universe existing  
	if (is.null(labels.universe)){ 
		mat.names <- sort(unique(unlist(dimnames(mat)))) 
	}else {mat.names<-labels.universe}
	
	
	#do we have all rows / columns?
	flag.row<-  mat.names %in% rownames(mat)
	flag.col<- mat.names %in% colnames(mat)
	
	mat.ext<-matrix(0, length(mat.names), length(mat.names) )
	rownames(mat.ext)<-colnames(mat.ext)<- mat.names	
	# fill extended matrix
	for (Row in rownames(mat)){
		for( Col in colnames(mat)){
			mat.ext[which(rownames(mat.ext) == Row),which(colnames(mat.ext) == Col)]<-mat[Row,Col]
		}
	}
return(mat.ext)
}
