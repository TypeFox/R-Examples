f.matrix.to.list <- function(x){
# CONVERTS A MATRIX x TO A LIST WHERE EACH COLUMN IS ONE ELEMENT OF THE
# LIST
#
	if(!(is.matrix(x) | is.data.frame(x))) stop("Object must be a matrix!")
	lapply(seq(dim(x)[2]), function(i, .mat) .mat[,i], .mat = x)

}
