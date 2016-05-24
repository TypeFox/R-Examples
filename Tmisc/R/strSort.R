#' Sort characters in a string
#' 
#' Alphabetically sorts characters in a string. Vectorized over x.
#' 
#' @author Stephen Turner
#' @keywords keywords
#'   
#' @param x A string to sort.
#'   
#' @return A sorted string.
#'   
#' @examples
#' strSort("cba")
#' strSort("zyxcCbB105.a")
#' strSort(c("cba", "zyx"))
#' strSort(c("cba", NA))
#' 
#' @export
strSort <- function(x) {
    ifelse(is.na(x), NA, sapply(lapply(strsplit(x, NULL), sort), paste, collapse=""))
}
