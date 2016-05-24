#  solve a generalized Cholesky matrix
solve.gchol <- function(a, b, full=TRUE, ...) {
    if (full) flag<-0 else flag<-1

    d <- a@Dim
    if (missing(b)) {
	# Return the inverse of the original matrix, for which a is the chol
	temp <- .C("gchol_inv", as.integer(d), 
                   x=as.double(a@.Data),
                   as.integer(flag))$x
	matrix(temp, d[1])
	}

    else {  # solve for right-hand side
	if (length(b) == d[1]) {
	    temp <- .C("gchol_solve", as.integer(d),
		                      as.double(a@.Data),
		                      y=as.double(b),
                                      as.integer(flag))$y
	    temp
	    }
	else {
	    if (!is.matrix(b) || nrow(b) != d[1]) 
		stop("number or rows of b must equal number of columns of a")
	    new <- b
	    for (i in 1:ncol(b)) {
		new[,i] <- .C("gchol_solve", as.integer(d),
			                     as.double(a@.Data),
		                             y=as.double(b[,i]),
                                             as.integer(flag))$y
		}
	    new
	    }
	}
    }
