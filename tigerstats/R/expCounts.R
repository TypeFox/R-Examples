#' @title Expected Cell Counts From Two-Way Tables

#' @description Simple instructional function to compute expected cell counts 
#' from a table of observed counts. 
#' 
#' 
#' @rdname expCounts
#' @usage expCounts(tab)
#' @param tab A table with two dimensions, or an object that can be coerced to one.
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @export
expCounts <- function(tab) {
  expected <- (rowSums(tab) %*% t(colSums(tab)))/sum(tab)
  rownames(expected) <- rownames(tab)
  colnames(expected) <- colnames(tab)
  return(expected)
}
