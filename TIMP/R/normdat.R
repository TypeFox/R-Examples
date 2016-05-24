"normdat" <-
function (mat) 
{
    if (!is.matrix(mat)) 
        mat <- as.matrix(mat)
    for (i in 1:ncol(mat)) 
    	if(max(abs(mat[, i])) == 0)
	  mat[, i] <- 1		 
        else
	  mat[, i] <- mat[, i]/max(abs(mat[, i]))
    mat
}
