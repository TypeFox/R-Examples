#' Validate n-gram
#' 
#' Checks if the character string may be used as an n-gram and its notation follows specific 
#' convention of \code{biogram} package.
#' @param x a \code{character} string representing single n-gram.
#' @return \code{TRUE} if n-gram's notation is correct, \code{FALSE} if not.
#' @export
#' @examples
#' print(is_ngram("1_1.1.1_0.0"))
#' print(is_ngram("not_ngram"))


is_ngram <- function(x) {
  sngram <- strsplit(x, "_")[[1]]
  
  if(!(length(sngram) %in% c(2, 3)))
    return(FALSE)
  
  #check if there is position information
  pos_inf <- ifelse(length(sngram) == 3, TRUE, FALSE)
  
  #validate position information
  if(pos_inf) 
    if(is.na(suppressWarnings(as.numeric(sngram[[1]]))) || as.numeric(sngram[[1]]) < 1)
      return(FALSE)
  
  seq <- strsplit(sngram[1 + pos_inf], ".", fixed = TRUE)[[1]]
  dists <- strsplit(sngram[2 + pos_inf], ".", fixed = TRUE)[[1]]
  
  #validate distance 
  if(length(seq) > 1) {
    #for n-grams biggers than uni-grams it must have length n - 1
    if(length(dists) != (length(seq) - 1)) 
      return(FALSE)
  } else {
    #for unigrams the distance vector has also size 1
    if(length(dists) != 1) 
      return(FALSE)
  }
  
  TRUE
}