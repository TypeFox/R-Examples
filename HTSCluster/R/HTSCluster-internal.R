.myloopfxn <- function(k, lambda, w.mat, s.mat, r, n, cols) {
	lambda.mat <- matrix(rep(rep(lambda[,k], times = r), each = n),
		nrow = n, ncol = cols)
	return(w.mat * s.mat * lambda.mat)
}
.myfxn <- function(var1, var2) {
	tmp <- var1 * log(var1/var2) + var2 - var1
	index <- which(var1 == 0)
	tmp[index] <- var2[index]
	rowSums(tmp)
}
.myprobafxn <- function(k, y, pi, mean) {
	pi[k] * exp(rowSums(dpois(y, mean[[k]], log=T)))
}
#.myrbind.fill.matrix <- function(matrices, rows, cols) 
#{
#	## Adapted from rbind.fill.matrix function in plyr package
#	## matrices is a list
#	## cols is the number of columns in each matrix
#	## rows is the number of rows in each matrix
#    cols <- 1:cols
#    rows <- rep(rows, length(matrices))
#    nrows <- sum(rows)
#    output <- matrix(NA, nrow = nrows, ncol = length(cols))
#    colnames(output) <- cols
#    pos <- matrix(c(cumsum(rows) - rows + 1, rows), ncol = 2)
#    for (i in seq_along(rows)) {
#        rng <- seq(pos[i, 1], length = pos[i, 2])
#        output[rng,] <- matrices[[i]]
#    }
#    output
#}
