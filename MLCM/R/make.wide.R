`make.wide` <-
function(d){
	nr <- nrow(d)
	nms <- substring(names(d)[1], 1, nchar((names(d)[1])) - 1)
	wts <- rep(c(1, -1), each = nr)
	ix.mat <- matrix(0, ncol = max(d), nrow = nr)
	ix.mat[matrix(c(rep(1:nr, 2), as.vector(unlist(d))), 
			ncol = 2)] <- wts
	ix.mat <- t(apply(ix.mat, 1, function(x) if (sum(x) == 0) x else
		rep(0, max(d))))
	colnames(ix.mat) <- paste(nms, seq(1, ncol(ix.mat)), sep = "")
	ix.mat[, -1, drop = FALSE]
}

`make.wide.full` <- function(d){
	cn <- function(y, mx) (y[, 2] - 1) * mxd[2] + y[, 1]
	nr <- nrow(d)
	mxd <- c(max(d[, 1:2]), max(d[, 3:4]))
	nc <- prod(mxd)	
#   Add column names
	nms <- sapply(seq(1, ncol(d), 2), function(x) 
		substring(names(d)[x], 1, 
		   nchar((names(d)[x])) - 1))
	fnm <- mapply(paste, nms, 
		list(seq_len(mxd[1]), seq_len(mxd[2])),
		sep = "", SIMPLIFY = FALSE)
	nms.f <- interaction(do.call(expand.grid, fnm), 
		sep = ":")
    ix.mat <- matrix(0, ncol = nc, nrow = nr)
    ix.mat[cbind(seq_len(nr), cn(d[, c(1, 3)], mxd))] <- 1
    ix.mat[cbind(seq_len(nr), cn(d[, c(2, 4)], mxd))] <- -1
    ix.mat <- t(apply(ix.mat, 1, function(x) if (sum(x) == 0) 
        x
    else rep(0, nc)))
    colnames(ix.mat) <- levels(nms.f)
	ix.mat[, -1, drop = FALSE]
	}
