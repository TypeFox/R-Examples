#' Count the total number of all common subsequences between a string and a vector/list of strings.
#' Subsequences cannot be interrupted by any item,
#' i.e. q-w is not considered a subsequence of q-e-w-r due to the interrupting 'e'
#'
#' @param vecA The single string
#' @param listB The vector/list of 1 or more strings
#' @param sep Delimiter separating each items in a sequence
#' @param dropFirstItem Boolean. If true, the first item in each sequence is excluded from counting all subsequences
#' @param ignoreLenOneSubseq Boolean. If true, all length one subequences are not counted as common subsequences
#' @param ignoreLenZeroSubseq Boolean. If true, the length zero subsequence (empty set) is not counted as a common subsequence
#' @return The total number of all common subsequences as an integer in a vector
#' @examples
#' calACSStrict("q-w-e-r", c("q-e-w-r","q-r-e-w"), "-")
#' calACSStrict("itemToBeDropped-q-w-e-r", "itemToBeDroped-q-e-w-r", "-", dropFirstItem=TRUE)
#'
#' @export

calACSStrict <- function(vecA, listB, sep="-", dropFirstItem=FALSE, ignoreLenOneSubseq=FALSE, ignoreLenZeroSubseq=FALSE){
  sapply(listB, calACSstrStrict, strA=vecA, sep=sep, dropFirstItem=dropFirstItem, ignoreLenOneSubseq=ignoreLenOneSubseq,
         ignoreLenZeroSubseq=ignoreLenZeroSubseq, simplify=TRUE, USE.NAMES = FALSE)
}
