#' Calculate the length of the longest common subsequence (KCS) between a string and a vector/list of strings.
#' Subsequences cannot be interrupted by any item,
#' i.e. q-w is not considered a subsequence of q-e-w-r due to the interrupting 'e'
#'
#' @param vecA The single string
#' @param listB The vector/list of 1 or more strings
#' @param sep Delimiter separating each items in a sequence
#' @param dropFirstItem Boolean. If true, the first item in each sequence is excluded from counting all subsequences
#' @return A list of vectors of the length of each common subsequence
#' @examples
#' lenACSStrict("q-w-e-r", c("q-e-w-r","q-r-e-w","q-w-r-e"), "-")
#' lenACSStrict("itemToBeDropped-q-w-e-r", "itemToBeDroped-q-e-w-r", "-", dropFirstItem=TRUE)
#'
#' @export

lenLCSStrict <- function(vecA, listB, sep="-", dropFirstItem=FALSE){
  unlist(sapply(listB, lenLCSstrStrict, strA=vecA, sep=sep, dropFirstItem=dropFirstItem, simplify = FALSE, USE.NAMES = FALSE))
}
