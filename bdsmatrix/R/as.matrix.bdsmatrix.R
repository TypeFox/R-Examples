as.matrix.bdsmatrix <- function(x, ...) {
    if (class(x) != 'bdsmatrix') stop('argument must be a bdsmatrix object')
    if (length(x@blocksize)==0) return(x@rmat)
    dd <- dim(x)
    d3 <- sum(x@blocksize)   # dim of square portion
    d4 <- sum(x@blocksize^2) # size of x@blocks
    newmat <- matrix(0., dd[1], dd[2], dimnames=x@Dimnames)
    temp <- .C('bdsmatrix_index1', 
	       as.integer(length(x@blocksize)),
	       as.integer(x@blocksize),
	       as.integer(c(1,0,0)),
	       as.integer(d3),
	       as.integer(1:d3 -1),
	       indexa = integer(d3*d3),
	       indexb = 0,
	       indexc = 0)$indexa

    newmat[1:d3, 1:d3] <- c(x@offdiag, x@blocks)[1+temp]
    if (length(x@rmat)>0) {
	newmat[, -(1:d3)] <- x@rmat
	newmat[-(1:d3),]  <- t(x@rmat)
	}
    newmat
    }

setAs('bdsmatrix', 'matrix', function(from)as.matrix.bdsmatrix(from))
setMethod('dim', 'bdsmatrix', function(x) x@Dim)
setMethod('dimnames', 'bdsmatrix', function(x) x@Dimnames)
setMethod('dimnames<-', 'bdsmatrix',
    function(x, value) {
	dd <- x@Dim
	if (is.null(value)) x@Dimnames <- NULL
	else {
	    if (is.list(value) && length(value)==2) {
		if (length(value[[1]])==0) val1 <- NULL
		else { 
		    val1 <- value[[1]]
		    if (length(val1) != dd[1]) 
			stop("Invalid length for row dimnames")
		    }
		if (length(value[[2]])==0) val2 <- NULL
		else { 
		    val2 <- value[[2]]
		    if (length(val2) != dd[2]) 
			stop("Invalid length for column dimnames")
		    }

		x@Dimnames <- list(val1, val2)
		}
	    else stop("dimnames must be a list of length 2")
	    }
	x
	})

print.bdsmatrix<- function(x, ...) print(as(x, 'matrix'), ...)
setMethod('show', 'bdsmatrix', 
	          function(object) show(as(object, 'matrix')))

setAs('bdsmatrix', 'vector', 
       function(from) as.vector(as.matrix.bdsmatrix(from)))

# this was commented out later: we don't want to inadvertently
#    create gigantic regular matrices 
#setIs('bdsmatrix', 'matrix', 
#      coerce=function(object)as.matrix.bdsmatrix(object))
      




