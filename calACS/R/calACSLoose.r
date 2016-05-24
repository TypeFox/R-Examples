#' Calculate the total number of all common subsequences between a string and a vector/list of strings.
#' Subsequences can be interrupted by items, i.e. q-w is considered a subsequence of q-e-w-r
#'
#' @param vecA The single string
#' @param listB The vector/list of 1 or more strings
#' @param sep Delimiter separating each items in a sequence
#' @param dropFirstItem Boolean. If true, the first item in each sequence is excluded from counting all subsequences
#' @return The total number of all common subsequences as an integer in a vector
#' @examples
#' calACSLoose("q-w-e-r", c("q-e-w-r","q-r-e-w"), "-")
#' calACSLoose("itemToBeDropped-q-w-e-r", "itemToBeDroped-q-e-w-r", "-", dropFirstItem=TRUE)
#'
#' @export

calACSLoose <- function(vecA, listB, sep="-", dropFirstItem=FALSE){
  sapply(listB, calACSstrLoose, strA=vecA, sep=sep, dropFirstItem=dropFirstItem, simplify=TRUE, USE.NAMES = FALSE)
}
