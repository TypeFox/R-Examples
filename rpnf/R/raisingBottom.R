#' returns true if given column c exceeds prevois column of same type (this is always column c-2)
#'
#' @param redData Data to consider
#' @param column Column to consider
raisingBottom <- function(redData,column) {
  # if (column-2>=1) # old version
  if (column>=3) # new, optimized version
    return (minBox(redData,column)>minBox(redData,column-2))
  else
    return (FALSE)
}
