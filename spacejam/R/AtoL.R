require("Matrix")
AtoL <-
function(A){
	p <- ncol(A)
	whA <- which(A)
	v1 <- whA%%p	
	v1[v1==0] <- p
	v2 <- ceiling(whA/p)
	ind <- v1 < v2
	
	return(as.matrix(rbind(as.integer(v1), as.integer(v2))[,ind]))
	}
