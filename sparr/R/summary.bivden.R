summary.bivden <- function(object, ...){

	print.bivden(x=object)
	
	if(length(object)==7) cat("Computed using ppp methods of package 'spatstat'\n")

	cat("Evaluated over",nrow(object$Zm),"by",ncol(object$Zm),"rectangular grid.\nDefined study region is a polygon with",length(vertices(object$WIN)$x),"vertices.\n\n")
	
	#cat("Estimated density description\n")
	
	cat("Estimated density range ",min(as.vector(object$Zm),na.rm=T)," to ",max(as.vector(object$Zm),na.rm=T),".\n",sep="")
	cat(sum(!is.na(as.vector(object$Zm))),"grid cells out of",prod(dim(object$Zm)),"fall inside study region.\n")
	
	#print(summary(as.vector(object$Zm)))

}
	