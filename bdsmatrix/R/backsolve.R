#
# The backsolve method for my matrices
#  If B= gchol(A) = LDL' the backsolve(B, x) solves L sqrt(D) y = x
# Since B is symmetric the transpose argument is ignored
#
# The next  lines are taken directly from the "Writing R Extensions" 
# manual. 
setGeneric("backsolve",
           function(r, ...) standardGeneric("backsolve"),
           useAsDefault= function(r, ...) base:::backsolve(r, ...))
           
#backsolve.default <- base:::backsolve
#formals(backsolve.default) <- c(formals(backsolve.default), alist(... = )) 

setMethod("backsolve", "gchol", 
    function(r, x, k = ncol(r), upper.tri=TRUE, ...) {
    if (any(diag(r) < 0))
        stop("Argument has a negative diagonal, cannot backsolve")

    if (!is.numeric(x)) stop("Invalid data type for x")
    x <- as.matrix(x)
    if (k!= floor(k)) stop("k must be an integer")
    if (k<1 || k > ncol(r)) stop("invalid value for k")

    if (nrow(x) != k)
        stop("Number of rows of x needs to match k")

    if (!is.logical(upper.tri) || is.na(upper.tri))
        stop("Invalid value for upper.tri option")
    storage.mode(x) <- "double"

    # I don't call with "r" itself, since the documentation on how
    #  to handle S4 classes internally is sparse to non-existent.
    # Looking at the code of Matrix, I can mimic, but don't trust.
    # The matrix x is fine though.
    drop(.Call("gcback", r@.Data, x, upper.tri, as.integer(k)))
})   
          
setMethod("backsolve", "gchol.bdsmatrix", 
    function(r, x, k=ncol(r), upper.tri=TRUE, ...) { 
          if (any(diag(r) < 0))
              stop("Argument has a negative diagonal, cannot backsolve")
          if (!is.numeric(x)) stop("Invalid data type for x")
          x <- as.matrix(x)
          if (k!= floor(k)) stop("k must be an integer")
          if (k<1 || k > ncol(r)) stop("invalid value for k")

          #Indexing a partial matrix would use less memory, but it's
          #  too much trouble in the remaining code.
          if (k < ncol(r)) r <- r[1:k, 1:k]

          if (nrow(x) != nrow(r))
              stop("Number of rows of x needs to match dimension of r")
          if (!is.logical(upper.tri) || is.na(upper.tri))
              stop("Invalid value for upper.tri optoin")
          storage.mode(x) <- "double"

          drop(.Call("gcback2", r@blocksize, r@blocks, r@rmat, 
                x, upper.tri))
      }) 


          
