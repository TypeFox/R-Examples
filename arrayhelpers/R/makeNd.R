##' Ensure/collapse an array into \code{n} dimensions and restore the old dimensions
##'
##' \code{nameNd} ensures a given number of dimensions:
##' If \code{a} has less than \code{N} dimensions, new dimensions of length 1 are appended.
##' If \code{a} has more than \code{N} dimensions, the supernumerary dimensions are collapsed onto
##' the last dimension.
##'
##' @param a an array (matrix, vector)
##' @param N the desired number of dimensions, 0 to remove the \code{dim} and \code{dimnames}
##' attributes (i.e. to create a vector). 
##' @return N-dimensional array
##' @author Claudia Beleites
##' @export  
##' @examples
##' v <- arrayhelpers:::v
##' v
##' makeNd (v, 1)
##' dim (makeNd (v, 1))
##' dim (makeNd (v, 3))
##' 
##' m <- arrayhelpers:::m
##' m
##' makeNd (m, 1)
##' dim (makeNd (m, 1))
##' makeNd (m, 0) 
##' dim (makeNd (m, 0))
##' makeNd (m, 3)
##' 
##' a <- arrayhelpers:::a
##' a
##' dim (makeNd (a, 1))
##' dim (makeNd (a, 0))
##' makeNd (a,  2)          
##' makeNd (a, -2)
##' makeNd (a, -4)
##' makeNd (a, 3);
##' 
makeNd <- function (a, N) {
   if (! all.equal (N, as.integer (N)))
     warning ("N truncated to integer value")
   N <- as.integer (N)

   d  <- dim (a)
   dn <- dimnames (a)

   ## push the old dimensions to the end of attribute old
   push (a, "old") <- list (list (names = attributes (a)$names, dimnames = dn, dim = d))
                                        # careful: names (a) and attr (a, "names") return dimnames
                                        # [[1]] for 1d arrays! --  We don't want that!
 
   if      (N == 0)        a <- .removedim   (a)
   else if (length (d) <  N && N > 0) a <- .appenddimafter    (a,  N, d, dn)
   else if (length (d) >  N && N > 0) a <- .collapsedimafter  (a,  N, d, dn)
   else if (length (d) < -N && N < 0) a <- .appenddimbefore   (a, -N, d, dn)
   else if (length (d) > -N && N < 0) a <- .collapsedimbefore (a, -N, d, dn)
    
   a
}

.test (makeNd) <- function (){
  checkEqualsNumeric (makeNd (v, 1), v)
  checkEqualsNumeric (makeNd (v, 2), v)
  checkEquals (dim (makeNd (v, 1)), 3L)
  checkEquals (dim (makeNd (v, 2)), c (3L, 1L))
  checkEquals (dimnames (makeNd (v, 1)) [[1]], names (v))

  checkEqualsNumeric (makeNd (m, 3), m)
  checkEquals (dim (makeNd (m, 3)), c(2L, 3L, 1L))
  checkEquals (dimnames (makeNd (m, 3)), c(dimnames (m), list (NULL)))
  
  checkEqualsNumeric (makeNd (a, 0), a)
  checkTrue (is.null (dim (makeNd (a, 0))))
  checkEquals (dim (makeNd (a, 2)), c (4L, 6L))
  checkEquals (dimnames (makeNd (a, 2)), list (rows = letters [1 : 4], NULL))

  checkTrue (is.null (attr (makeNd (ensuredim (v), 2), "old")[[1]]$names),
             msg = "names attribute for 1d arrays")
}

.removedim <- function (a){
  dim (a)      <- NULL
  dimnames (a) <- NULL

  a
}

.appenddimafter <- function (a, N, d, dn){
  if (is.null (d)){   # vector
    d <- length (a)
    dn <- list (names (a))
  }
   
  dim (a)        <- c (d,  rep (1,           N - length (d)))
  if (! is.null (dn))
    dimnames (a) <- c (dn, rep (list (NULL), N - length (dn)))
  
  a
}
.appenddimbefore <- function (a, N, d, dn){
  if (is.null (d)){   # vector
    d <- length (a)
    dn <- list (names (a))
  }
   
  dim (a)        <- c (rep (1,           N - length (d )),  d)
  if (! is.null (dn))
    dimnames (a) <- c (rep (list (NULL), N - length (dn)), dn)
  
  a
}

.collapsedimafter <- function (a, N, d, dn) {
  dim (a)        <- c (d  [seq_len (N - 1)], prod (d [N : length (d)]))
  if (! is.null (dn))
    dimnames (a) <- c (dn [seq_len (N - 1)], list (NULL))

  a  
}

.collapsedimbefore <- function (a, N, d, dn) {
  dim (a)        <- c (prod (d [seq_len (length (d) - N + 1)]), tail (d,  N - 1))
  if (! is.null (dn))
    dimnames (a) <- c (list (NULL),                             tail (dn, N - 1))

  a  
}
