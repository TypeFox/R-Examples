#' Interleave
#'
#' Used Internally
#' Interleave two character vectors. The output is a vector. The first entry is from the first vector.
#' Vectors can be of different lengths. If one is shorter than the other, entries of unmatched longer vector are left un-interleaved.  
#' @param vec1 first vector
#' @param vec2 second vector
#' @return interleaved vector
#' # t1 <- paste0("t1", letters[1:5]); t2 <- paste0("t2", letters[1:5]); interleave(t1, t2)

interleave <- function(vec1, vec2) 
{
  if (!is.vector(vec1, mode="character") | !is.vector(vec2, mode="character")) stop("Either or both of the inputs are not a character vector.\n The function only takes character vectors.")
  if (length(vec1) != length(vec2)) cat("Vectors are of different lengths. Unmatched entries of the longer vector will be left un-interleaved.\n")
  
  ord1 <- 2*(1:length(vec1))-1
  ord2 <- 2*(1:length(vec2))
  
  c(vec1,vec2)[order(c(ord1,ord2))]
}
