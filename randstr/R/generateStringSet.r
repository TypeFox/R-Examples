generateStringSet <- function(nString, len, symbolSet, delimiter){
  stringSet <- ""
  for(iString in 1:nString){
    curString <- ""
    for(iItem in 1:len[iString]){
      symboln <- sample((1:length(symbolSet)), 1)
      symbol <- symbolSet[symboln]
      curString <- paste(curString, symbol, sep = delimiter)
    }
    curString <- stri_sub(curString, from = length(delimiter) + 1)
    stringSet[iString] <- curString
  }
  stringSet
}