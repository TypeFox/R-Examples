##' \code{restoredim} restores the shape the array had before \code{makeNd} was called.
##'
##' Attributes \code{old.dim} and \code{old.dimnames} are used by default. \code{restoredim} is the
##' inverse of \code{makeNd}.
##'
##' Note that missing attributes as well as \code{old.dim = NULL} produce a (dimensionless)
##' vector. This is also the case if \code{a} lost the \code{old.*} attributes during 
##' computations like \code{as.numeric}, \code{c}, etc..
##'
##' \code{fromend} together with numeric \code{usedim} specifies dimensions counting from the
##' end. E.g. \code{fromend = TRUE} and \code{usedim = 1 : 3} for an array to be restored to 10d
##' means restoring dimensions 8 : 10. \code{fromend = TRUE} and \code{usedim = -(1 : 3)} restores
##' dimensions 1 to 7.
##' @param old list containing a list with (possibly) elements \code{dim}, \code{dimnames}, and
##' \code{names}. The nth last element of this list is used.
##' @param n how many makeNdim steps to go back?
##' @param ... ignored
##' @param usedim use only the specified dimensions
##' @param fromend if \code{TRUE}, numeric \code{usedim} are counted from the end, see details.
##' @param drop should 1d arrays drop to vectors?
##' @return an array
##' @author Claudia Beleites
##' @rdname makeNd
##' @export
##' @examples
##'
##' a <- array (1 : 24, 4 : 3)
##' a
##' restoredim (makeNd (a, 0))
##'
##' x <- makeNd (a, 0)
##' attr (x, "old")
##' 
restoredim <- function (a, old = NULL, n = 1L, ...,
                        usedim = TRUE, fromend = FALSE, drop = FALSE){

  if (is.null (old)){
    old <- peek (a, "old", n = n)
    a <- pop (a, "old", n = n)
  } else if (!is.null (names (old [[n]])) &&  all (names (old [[n]]) %in% c("dim", "dimnames", "names")) &&
             (is.null (names (old      )) || !all (names (old      ) %in% c("dim", "dimnames", "names")))
             )
      old <- old [[n]]
  ## else assume a list with dim, dimnames and names
  
  if (fromend && is.numeric (usedim))
    usedim <- sort (rev (seq_along (old$dim)) [usedim])
  else
    usedim <- sort (numericindex (x = old$dim, i = usedim, n = names (old$dimnames)))

  a <- structure (a, .Dim = old$dim [usedim],
                  .Dimnames =  lon (old$dimnames [usedim]),
                  .Names =  old$names)

  drop1d (a, drop = drop)
}

.test (restoredim) <- function (){
  checkIdentical (a, restoredim (makeNd (a,  5)))
  checkIdentical (a, restoredim (makeNd (a,  0)))
  
  checkIdentical (a, restoredim (makeNd (a, -5)))
  checkIdentical (a, restoredim (makeNd (a, -2)))

  checkIdentical (v, restoredim (makeNd (v,  0)))
  checkIdentical (v, restoredim (makeNd (v,  1)))
  checkIdentical (v, restoredim (makeNd (v, -1)))
  
  checkIdentical (dim (restoredim (as.numeric (makeNd (a, 0)))), NULL)

  checkIdentical (a, restoredim (restoredim (makeNd (makeNd (a,  5), 0))))
  checkIdentical (a, restoredim (makeNd (makeNd (a,  5), 0), n = 2))

  checkIdentical (a, restoredim (makeNd (makeNd (a,  5), 0), n = 3)) # OK

  warn <- options(warn = 2)$warn
  on.exit (options (warn = warn))
  checkException (restoredim (makeNd (makeNd (a,  5), 0), n = 3))

  tmp <- makeNd (v, 5)
  old <- attr (tmp, "old")
  tmp <- pop (tmp, "old")
  checkIdentical (v, restoredim (tmp, old = old [[1]]))
  checkIdentical (v, restoredim (tmp, old = old))

  ## TODO: test drop and usedim
}

