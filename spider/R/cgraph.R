cgraph <- function(x, y = NULL, ...){
	if(!is.null(y)) mat <- cbind(x, y) else mat <- x
	dd <- dim(mat)[1]
	if(dd < 2) return()
	if(is.null(dd)) return()
	ddComb <- combn(1:dd, 2)
	segments(mat[ddComb[1,], 1], mat[ddComb[1,], 2], mat[ddComb[2,], 1], mat[ddComb[2,], 2], ...)
}
