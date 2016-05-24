generateRandSymbolSet <- function(symbolSetSize){
  if(symbolSetSize <= 26){
    symbolSet <- randomStrings(n = symbolSetSize, len = 1, digits = FALSE, upperalpha = FALSE, loweralpha = TRUE, 
                               unique = TRUE)
  }
  else if(symbolSetSize <= 52){
    symbolSet <- randomStrings(n = symbolSetSize, len = 1, digits = FALSE, upperalpha = TRUE, loweralpha = TRUE, 
                               unique = TRUE)
  }
  else if(symbolSetSize <= 62){
    symbolSet <- randomStrings(n = symbolSetSize, len = 1, digits = TRUE, upperalpha = TRUE, loweralpha = TRUE, 
                               unique = TRUE)
  }
  else{
    symbolSet <- "err in generateRandSymbolSet"
  }
}
