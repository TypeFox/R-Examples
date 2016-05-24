#' IBD probabilities
#' 
#' Returns the IBD probabilities for a pairwise relationship.
#' 
#' @param x A pairwise relationship (e.g. "FS" for full siblings), or a 3-vector of IDB-probabilities.
#' @details When a character vector is supplied, the function returns the IBD probabilities. The following pairwise relationships are supported:
#' \itemize{
#' \item "ID", identity/monozygotic twins
#' \item "FS", full siblings
#' \item "PO" alias "PC", parent/offspring
#' \item "HS", half siblings
#' \item "AV", uncle/nephew
#' \item "FC", first cousins
#' \item "SC", second cousins
#' \item "UN", unrelated
#' }
#' 
#'          When a numeric vector is supplied (e.g. c(0,1,0)), the function checks the input for validity and returns the input; throws an error otherwise.
#' @return Numeric of length 3 with the probabilities that 0, 1 or 2 alleles are identical by descent
#' @examples identical(ibdprobs("PO"),ibdprobs(c(0,1,0))) #TRUE
#' @export
ibdprobs <- function(x){
  if (missing(x)) stop("Please supply IBD-probabilities either as a length 3 numeric or a character vector, e.g. \"FS\" or \"UN\"")
    
  if (is.numeric(x)){
    if (any(length(x)!=3L,!all.equal(sum(x),1.),x<0,x>1)) stop("IBD probabilities should be a numeric vector of length 3 with sum 1")
  }else if (is.character(x)){
    if (length(x)==1){
      if (is.null(Zibdpr[[x]])) stop("Unknown hypothesis: ", x,". Choose one of ",paste(names(Zibdpr),collapse=", "),".")      
      return(ibdprobs(Zibdpr[[x]]))
    }else{
      stop("x should have length 1")
    } 
    
  }else{
    stop("x should be either a numeric vector of IBD probabilities or a character vector indicating the hypothesis")
  }
  x
}