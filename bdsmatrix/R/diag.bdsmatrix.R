setMethod('diag', 'bdsmatrix', function(x, nrow, ncol) {
    if (class(x) != 'bdsmatrix') stop("Argument must be a bdsmatrix object")
    
    d <- x@Dim
    d3 <- sum(x@blocksize)
    temp <- .C('bdsmatrix_index1',
	       as.integer(length(x@blocksize)),
	       as.integer(x@blocksize),
	       as.integer(c(0,1,0)),
	       as.integer(d3),
	       as.integer(1:d3 -1),
	       integer(1),
	       indexb = integer(d3),
	       integer(1))$indexb

    if (length(x@rmat) > 0) {
	temp2 <- seq(from=d3+1, by= d[2]+1, length= d[1] - d3)
	c(x@blocks[temp], x@rmat[temp2])
	}
    else x@blocks[temp]
    })

setMethod("diag<-","bdsmatrix" ,function(x, value) {
    if (class(x) != 'bdsmatrix') stop("Argument must be a bdsmatrix object")
    
    d <- x@Dim
    if (length(value) != d[1]) stop("Wrong length for diagonal")
    d3 <- sum(x@blocksize)
    temp <- .C('bdsmatrix_index1',
	       as.integer(length(x@blocksize)),
	       as.integer(x@blocksize),
	       as.integer(c(0,1,0)),
	       as.integer(d3),
	       as.integer(1:d3 -1),
	       integer(1),
	       indexb = integer(d3),
	       integer(1))$indexb
    x@blocks[temp] <- value[1:d3]
    if (length(x@rmat) > 0) {
	temp2 <- seq(from=d3+1, by= d[2]+1, length= d[1] - d3)
	x@rmat[temp2] <- value[-(1:d3)]
	}
    x
    })
    
