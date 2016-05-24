#' Pretty Printing
#' 
#' Prettier printing for matrices and data frames.
#' 
#' @param x An object of class \code{"matrix"} or \code{"data.frame"}.
#' @param rowdots Integer specifying the row to replace with \code{...} 
#'   notation. Default is 4.
#' @param coldots Integer specifying the column to replace with \code{...} 
#'   notation. Default is 4.
#' @param digits The minimum number of significant digits to be printed in 
#'   values.
#' @param ... Additional optional arguments. None are used at present.
#' 
#' @details 
#' For object of class \code{"matrix"} or \code{"data.frame"} (which are coerced 
#' to a matrix via the \code{data.matrix} function), \code{pprint} will replace 
#' all the rows starting from \code{rowdots} up to and including the second-to-last 
#' row with a single row filled with \code{...}s. The same is applied to the 
#' columns as well. Hence a large matrix (or data frame) will be printed in a 
#' much more compact form. 
#' @export
#' @examples
#' pprint(randn(100, 100))
#' pprint(resize(1:100, 10, 10))
pprint <- function(x, ...) {
  UseMethod("pprint")
}


#' @rdname pprint
#' @method pprint matrix
#' @export
pprint.matrix <- function(x, rowdots = NULL, coldots = NULL, digits = NULL, 
                          ...) {
  
  # Default values
  if (is.null(rowdots)) rowdots <- getOption("pprint.rowdots")
  if (is.null(coldots)) coldots <- getOption("pprint.coldots")
  if (is.null(digits)) digits <- getOption("digits")
  
  # Row labels
  row_labels <- if (is.null(rownames(x))) {
    paste0("[", seq_len(nrow(x)), ",]")
  } else {
    rownames(x)
  }
  
  # Columns labels
  col_labels <- if (is.null(colnames(x))) {
    paste0("[,", seq_len(ncol(x)), "]")
  } else {
    colnames(x)
  }
  
  # Convert to character matrix (after rounding, if appropriate)
  charx <- if (typeof(x) == "character") {
    x
  } else if (typeof(x) %in% c("integer", "logical")) {
    as.character(x)
  } else {
    sprintf(paste0("%.", digits, "f"), x)
  }
  dim(charx) <- dim(x)
  
  # Case 1: rows and columns do not have dots
  if (nrow(x) <= rowdots + 1 && ncol(x) <= coldots + 1) {
    res <- x  
  }
  
  # Case 2: rows have dots, columns do not
  if (nrow(x) > rowdots + 1 && ncol(x) <= coldots + 1) {
    res <- rbind(as.matrix(charx[seq_len(rowdots - 1), ]), 
                 rep("...", ncol(charx)),
                 charx[nrow(charx), ])
    row_labels <- add_dots(row_labels, pos = rowdots) 
  }
  
  # Case 3: rows do not have dots, columns have dots
  if (nrow(x) <= rowdots + 1 && ncol(x) > coldots + 1) {
    res <- t(apply(charx, 1, add_dots, pos = coldots))
    col_labels <- add_dots(col_labels, pos = coldots)
  }
  
  # Case 4: rows and columns have dots
  if (nrow(x) > rowdots + 1 && ncol(x) > coldots + 1) {
    # Add first rowdots-1 rows
    smallx <- t(apply(charx[seq_len(rowdots - 1), ], 1, add_dots, 
                      pos = coldots))
    res <- rbind(smallx, 
                 rep("...", ncol(smallx)),
                 add_dots(charx[nrow(charx), ], pos = coldots))
    row_labels <- add_dots(row_labels, pos = rowdots)
    col_labels <- add_dots(col_labels, pos = coldots)
  } 
  
  # Print "pretty" matrix
  cat(desc_mat(x), "\n")
  cat("\n")
  prmatrix(res, rowlab = row_labels, collab = col_labels, quote = FALSE, 
           right = TRUE)
  
  # Return a (temporarily) invisible copy of x
  invisible(x)
  
}


#' @rdname pprint
#' @method pprint data.frame
#' @export
pprint.data.frame <- function(x, rowdots = NULL, coldots = NULL, digits = NULL, 
                              ...) {
  pprint(data.matrix(x), rowdots = rowdots, coldots = coldots, digits = digits,
         ...) 
}