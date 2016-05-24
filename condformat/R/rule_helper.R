
#' Returns a matrix with the value of the CSS field for each element
#' of our data frame.
#'
#' @param finalformat list with the css_fields
#' @param field the name of the CSS field to be returned
#' @return a matrix of size xview with the CSS field information
get_css_field <- function(finalformat, field) {
  if (!(field %in% names(finalformat$css_fields))) {
    field_matrix <- matrix(data = "",
                           nrow = nrow(finalformat$css_cell),
                           ncol = ncol(finalformat$css_cell))
  } else {
    field_matrix <- finalformat$css_fields[[field]]
  }
  return(field_matrix)
}

#' If rule_locks is TRUE, locks the cells specified by mask
#' so no further rules modify them
#'
#' @param rule_locks a logical. If FALSE, lock_cells does nothing
#' @param mask a logical matrix. TRUE if the cell is affected by the rule
#' @param finalformat The finalformat list, with the format options
#'
#' @return finalformat, the list with format options
#'
lock_cells <- function(rule_locks, mask, finalformat) {
  if (identical(rule_locks, TRUE)) {
    finalformat$css_cell_unlocked[mask] <- FALSE
  }
  return(finalformat)
}

fill_css_field_by_cols <- function(finalformat, field, values, columns, xview, lockcells) {
  index.j <- match(columns, colnames(xview))
  index.j <- index.j[!is.na(index.j)]

  mask <- finalformat$css_cell_unlocked
  # mask == TRUE if cell can be changed, false otherwise
  # We can't change columns not affected:
  mask[,-index.j] <- FALSE

  backgr <- get_css_field(finalformat, field)
  backgr[mask] <- values[mask]

  finalformat$css_fields[[field]] <- backgr
  finalformat <- lock_cells(lockcells, mask, finalformat)
  return(finalformat)
}

