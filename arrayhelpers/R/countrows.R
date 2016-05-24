
##' Count equal rows
##'
##' matrices are converted to data.frame.
##' @param x the matrix or data.frame
##' @return data frame with unique rows, their counts and indices into the original data.frame
##' @note this function is subject to changes in the future.
##' @author Claudia Beleites
##' @rdname countrows
##' 
countRows <- function(x) {
## TODO: derive matrix function
## TODO: generalize to nd arrays
## TODO: tests
## TODO: examples
  
  if (is.matrix (x) && (ncol (x) == 1))
    x <- as.vector (x) 

  order.x <- do.call (order, as.data.frame (x))
  
  if (is.vector (x)) {
    equals.prev <-          x [tail (order.x, -1)  ] == x [head (order.x, -1)  ]
  } else {
    equals.prev <- rowSums (x [tail (order.x, -1), ] != x [head (order.x, -1), ]) == 0
  }

  indices <- split (order.x, cumsum (c (TRUE, ! equals.prev)))

  if (is.vector (x)) {
    x <- x [sapply (indices, function (v) v [[1]]),   drop = FALSE]
  } else {
    x <- x [sapply (indices, function (v) v [[1]]), , drop = FALSE]
  }

  data.frame (cts = sapply (indices, length),
              ind = I (indices),
              x)
}
