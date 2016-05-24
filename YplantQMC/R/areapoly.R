areapoly <- function(x.mat){
	  x.segmat <- cbind(x.mat, rbind(x.mat[2:nrow(x.mat), ],x.mat[1, ]))
    abs(sum(x.segmat[,1] * x.segmat[,4] - x.segmat[,3] * x.segmat[,2])) / 2
}		