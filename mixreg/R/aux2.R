aux2 <- function(t11,t12,t13,t21,t22,t23,t31,t32,t33) {
#
# Auxilliary function to pack the bits and pieces of
# info2 into a matrix.
#
	rslt <- cbind(rbind(t11,t21),c(t12,t22))
	if(!is.null(t13)) {
		rslt <- cbind(rslt,c(t13,t23))
		if(is.null(t31)) return(rslt)
		else return(rbind(rslt,c(t31,t32,t33)))
	}
	if(is.null(t31)) return(rslt)
	rbind(rslt,c(t31,t32))
}
