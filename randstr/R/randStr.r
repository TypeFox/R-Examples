#' Generate random strings of a dictated size of symbol set and distribution of the lengths of strings.
#' The output is a list of 2 items: [[1]]vector of strings, [[2]]vector of the corresponding symbol set
#'
#' @import truncnorm
#' @import random
#' @import stringi
#'
#' @param nString Number of strings to be generated
#' @param maxLen Maximum allowable length of the strings
#' @param minLen Minimum allowable length of the strings
#' @param stdevLen Standard deviation of the length of the strings, applicable if using truncated normal distribution for the lengths of strings i.e. distLen='normal'
#' @param distLen Distribution of the lengths of strings, can take on value of 'normal' or 'uniform'
#' @param symbolSetSize Allowable number of different symbols to appear in the strings
#' @param delimiter symbol separating each item in the strings
#'
#' @examples randStr(nString = 10, maxLen = 30, minLen = 1, stdevLen = 15, symbolSetSize = 25)
#' @examples randStr(nString = 10, maxLen = 30, minLen = 1, symbolSetSize = 25, distLen = 'uniform')
#' @export

randStr <- function(nString, maxLen, minLen, stdevLen = 0, distLen = 'normal', symbolSetSize, delimiter = "-"){
#   require(truncnorm)
#   require(random)
#   require(stringi)

  if(!nString%%1==0 || !maxLen%%1==0 || !minLen%%1==0 || !symbolSetSize%%1==0){
    return("err: nString, maxLen, minLen or symbolSetSize not integer")
  }
  if(nString<0 || maxLen<0 || minLen<0 || symbolSetSize<0 || stdevLen<0){
    return("err: nString, maxLen, minLen, symbolSetSize or stdevLen is negative")
  }
  if(distLen != 'normal' && distLen != 'uniform'){
    return("err: typo in distLen")
  }
  if(symbolSetSize > 62){
    return("err: symbol set size bigger than the number of characters available (lower and upper case alphabets and numbers 0-9)")
  }

  len <- generatelens(maxLen, minLen, stdevLen, nString, distLen)
  symbolSet <- generateRandSymbolSet(symbolSetSize)
  stringSet <- generateStringSet(nString, len, symbolSet, delimiter)

  list(stringSet, symbolSet)
}

randStr(nString = 10, maxLen = 30, minLen = 1, stdevLen = 15, symbolSetSize = 25, delimiter = "-", distLen = 'normal')
