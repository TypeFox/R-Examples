#
# Code for the generalized cholesky  A = LDL', where L is lower triangular
#   with 1's on the diagonal, and D is diagonal.
# The decompostions exists for any square symmetric matrix.  
# If A is positive definite, then all elements of D will be positve.
# If A is not full rank, then 0's on the diagonal of D signal the redundant 
#   columns.  
# Note that gchol is both a class (setClass) and a generic function.
#
setClass('gchol', 
	 representation(.Data= 'numeric',
			Dim = 'integer',
			Dimnames = 'list',
			rank = 'integer'))
			
setGeneric('gchol', function(x, tolerance=1e-10) standardGeneric('gchol'),
           useAsDefault=FALSE)

as.matrix.gchol <- function(x, ones=TRUE, ...) {
    temp <- matrix(x@.Data, x@Dim[1], dimnames=x@Dimnames, byrow=TRUE)
    if (ones) diag(temp) <- 1
    temp
    }

setAs('gchol', 'matrix', function(from) as.matrix.gchol(from))

setMethod('gchol', signature(x='matrix'),
    function(x,  tolerance) {
	d <- dim(x)
	if (d[1] != d[2]) 
		stop("Cholesky decomposition requires a square matrix")
#	if (!is.logical(all.equal(as.vector(x), as.vector(t(x)))))
#		stop("Cholesky decomposition requires a symmetric matrix")
	temp <- .C("gchol", as.integer(d[1]),
		   x =   as.double(x),
		   rank= as.double(tolerance))
      
        dnames <- dimnames(x)
        if (is.null(dnames)) dnames <- list(NULL, NULL)

	new('gchol', .Data= temp$x  , Dim=d, 
	    Dimnames= dnames, rank=as.integer(temp$rank))
	})

setMethod('diag', signature(x='gchol'),
    function(x, nrow, ncol) {
	d <- x@Dim[1]
	x@.Data[ seq(1, length=d, by=d+1)]
	})


setMethod('show', 'gchol', function(object) show(as.matrix(object, F)))
setMethod('dim', 'gchol',  function(x) x@Dim)
setMethod('dimnames', 'gchol', function(x) x@Dimnames)
 	

# Multiplication methods.
# If the gchol can be written as Cholesky decompostion, i.e., if
#   all of the diagonal elements are >=0, then return the product
#   of the cholesky with the vector or matrix.  Otherwise squawk.
# 
setMethod("%*%", signature(x='gchol', y='matrix'),
    function(x, y) {
        if (!is.numeric(y))
            stop("Matrix multiplication is defined only for numeric objects")
        dy <- dim(y)
        dx <- dim(x)
        ldy <- length(dy)
        if (ldy!=2) dy <- c(length(y), 1)
        if (dx[2] != dy[1]) 
            stop("Number of columns of x should be the same as number of rows of y")
        
        if (any(diag(x) < 0)) stop("gchol matrix does not have a Cholesky repres
entation: no matrix product is possible")

        as.matrix(x) %*% (sqrt(diag(x)) * y)
        })

setMethod("%*%", signature(x='matrix', y='gchol'),
    function(x, y) {
        if (!is.numeric(x))
            stop("Matrix multiplication is defined only for numeric objects")
        dy <- dim(y)
        dx <- dim(x)
        ldx <- length(dx)
        if (ldx!=2) dx <- c(length(x), 1)
        if (dx[2] != dy[1]) 
            stop("Number of columns of x should be the same as number of rows of y")
        
        if (any(diag(y) < 0)) stop("gchol matrix does not have a Cholesky repres
entation: no matrix product is possible")

        (y %*% as.matrix(x)) * rep(sqrt(diag(y)), each=ncol(y))
        })

setMethod('[', "gchol", 
 function(x, i,j, drop=TRUE) {
     if (missing(i) && missing(j)) return(x)
     temp <- matrix(x@.Data, nrow=x@Dim[1], dimnames=x@Dimnames)
     if (missing(i)) temp[,j,drop=drop]
     else {
         if (missing(j)) temp[i,,drop=drop]
         else {
             temp <- temp[i,j,drop=drop]
             if (length(i)==length(j) && length(i)>1 && all(i==j)) {
                 # in this case only, the result is a gchol object
                 new("gchol", .Data= as.vector(temp), Dim=dim(temp),
                     Dimnames=dimnames(temp), rank=sum(diag(temp) !=0))
             }
             else temp
         }
     }
})


 
