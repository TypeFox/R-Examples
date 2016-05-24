#' Matrices
#'
#' Like \code{matrix}, \code{mat} creates a matrix from the given set of 
#' values. However, these values can also be represented by a character string, 
#' or a list of vectors.
#' 
#' @param x A data vector, character string, or a list.
#' @param rows Logical. If TRUE (the default) the matrix is filled by rows, 
#'             otherwise the matrix is filled by columns.
#' @param sep Separator string. Values within each row/column of x are 
#'            separated by this string. Default is \code{","}.
#' @param ... Aditional optional arguments.
#' @return A matrix of class \code{c("matrix", "mat")}.
#' @seealso \code{\link{bmat}}, \code{\link{dmat}}, \code{\link{matrix}}.
#' @export
#' @examples
#' ## Using character vectors
#' mat("1, 2, 3, 4; 5, 6, 7, 8")  # ";" separates rows
#' mat("1, 2, 3, 4; 5, 6, 7, 8", rows = FALSE)  # ";" separates columns
#' mat("1 2 3 4; 5 6 7 8", sep = "")  # use spaces instead of commas
#' mat(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2, byrow = TRUE)  # works like matrix too
#'
#' ## Using a list
#' z1 <- list(1:5, 6:10)
#' z2 <- list(a = 1:5, b = 6:10)
#' mat(z1)
#' mat(z2)  # preserves names as row names
#' mat(z2, rows = FALSE)  # preserves names as column names
mat <- function(x, ...) {
  UseMethod("mat")
}


#' @rdname mat
#' @method mat default
#' @export
mat.default <- function(x, ...) {
  matrix(x, ...)  # default to base matrix function
}


#' @rdname mat
#' @method mat character
#' @export
mat.character <- function(x, rows = TRUE, sep = getOption("mat.sep"), 
                          ...) {
  
  ## Gather rows and individual values
  # seps <- paste0("[;", "sep", "]")  # separate all at once
  vecs <- unlist(strsplit(x, split = ";"))  # column/row vectors
  char_vals <- if (!is.null(sep)) {
    unname(unlist(lapply(vecs, strsplit, split = sep)))
  } else {
    vecs
  }
  num_vals <- unlist(lapply(char_vals, function(x) eval(parse(text = x))))
  
  ## Form matrix from parsed values by calling R's built-in matrix function
  if (rows) {
    matrix(num_vals, nrow = length(vecs), byrow = TRUE, ...)
  } else {
    matrix(num_vals, ncol = length(vecs), byrow = FALSE, ...)
  }
  
}


#' @rdname mat
#' @method mat list
#' @export
mat.list <- function(x, rows = TRUE, ...) {
  
  ## Check element types
  if (!all(sapply(x, class) %in% c("numeric", "integer"))) {
    stop("Each element must be of type 'numeric' or 'integer'.", call. = FALSE)
  }
  
  ## Check length of each element
  if (!all(sapply(x, length) >= 1) && length(unique(sapply(x, length))) != 1) {
    stop("Each element must contain at least one value.", call. = FALSE)
  }
  
  ## Form matrix by combining elements
  if (rows) do.call(rbind, x) else do.call(cbind, x)
  
}
