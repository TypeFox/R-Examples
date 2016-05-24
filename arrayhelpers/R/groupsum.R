##' Extension of \code{rowsum}
##' 
##' \code{groupsum} extends \code{\link{rowsum}}: it allows \code{group} to be an array of the same shape
##' as \code{x}. 
##' @param x array to be \code{rowsum}med
##' @param group grouping variable (integer or factor) indicating groups of samples. \code{}
##' @param dim along which dimension should the group sums be taken? (default: rows)
##' @param reorder should the groups be ordered? see \code{\link[base]{rowsum}}
##' @param na.rm shoud \code{NA}s be removed?
##' @param ... ignored
##' @param drop should 1d arrays drop to vectors?
##' @return like \code{\link[base]{rowsum}}, but further dimensions of the array are preserved.
##' @author Claudia Beleites
##' @seealso \code{\link[base]{rowsum}} \code{\link{rowsum}}
##' @keywords array algebra arith
##' @export
groupsum <- function(x, group = NULL, dim = 1L, reorder=TRUE, na.rm = FALSE, ...,
                     drop = ! is.array (x)) {
  x <- ensuredim (x)

  ## permute the group dimension to the beginning
  x <- aperm (x, c (dim, seq_len (ndim (x)) [-dim]))

  x <- makeNd (x, 2)
  old <- attr (x, "old") [[1]]
  x <- pop (x, "old")
  

  if (is.null (group)){                 # almost no gain...
    x   <- colSums (x, na.rm = TRUE, drop = FALSE)
  } else if (is.null (dim (group))) {
    if (length (group) == nrow (x)) {
      x   <- rowsum  (x, group = group, na.rm = TRUE)
    } else {
      stop ("wrong length: group ", length (group), " dim: ", nrow (x))
    }   
  } else if (all (dim (group) == dim (x))) {
    stop ("grouping factors of size dim (p) not yet implemented.")
  } 

  old$dim [1] <- nrow (x)
  old$dimnames [1] <- list (rownames (x))

 # mostattributes (x) <- old

  x <- restoredim (x, old = old)
  x <- aperm (x, order (c (dim, seq_len (ndim (x)) [-dim])))

  drop1d (x, drop = drop)
}

.test (groupsum) <- function (){
  groups <- c(2, 1, 2)
  checkEquals (groupsum (a, group = groups, dim = 2),
               structure(c(5 : 8, (5 : 8) * 2L, 17 : 20, (17 : 20) * 2L), 
                           .Dim = c(4L, 2L, 2L),
                           .Dimnames = structure(list(rows = letters [1:4],
                             columns = c("1", "2"),
                             d3 = c("1", "2")),
                             .Names = c("rows", "columns", "d3"))))

  b <- a
  dim (b) <- c (2, 2, 3, 2)
  g <- groupsum (b, group = groups, dim = 3)

  checkEquals (b [,,1,] + b [,,3,], g [,,2,])
  checkEquals (b [,,2,]           , g [,,1,])
  
}
