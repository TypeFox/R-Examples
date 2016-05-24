convUserFunctions <- function(linesMat){
  funInd <- which(grepl("^function", linesMat))
  linesFun <- linesMat[funInd]

  allEndInd <- which(grepl("^}", linesMat))
  retInd <- allEndInd[vapply(funInd[-1], function(x){
    rev(which(allEndInd < x))[1]
  }, 1)] -1
  retInd <- c(retInd, rev(allEndInd)[1] - 1 )

  #Deal with return objects
  returnObj <- getBetween(linesFun, "[", "]")
  multOutSet <- grepl(',', returnObj)

  returnObj[multOutSet] <- sprintf("list(%s)", returnObj[multOutSet])
  returnObj <- sprintf("\treturn(%s)", returnObj)

  #Deal with the arguments
  args <- getBetween(linesFun, "(", ")")
  varargSet <- grepl("varargin", args)
  comSet <- grepl("^#", linesMat) | grepl("^\\s+#", linesMat)
  noComInd <- which(!comSet)
  decInd <- noComInd[vapply(funInd, function(x){
    which(noComInd > x)[1]
    }, 1)] - 1
  decObj <- rep("\tvarargin <- list(...)", length(which(varargSet)))
  decInd <- decInd[varargSet]


  funName <- getBetween(linesFun, "- ", "(")
  linesFun <- paste(funName, "<-", sprintf("function(%s){", args))
  linesFun <- gsub("varargin", "...", linesFun)

  linesMat[funInd] <- linesFun

  # save the appending stuff for last cause it throws off all the indices
  allAppObj <- c(returnObj, decObj)
  allAppInd <- c(retInd, decInd)

  for (funi in rev(order(allAppInd))){
  	linesMat <- append(linesMat, allAppObj[funi],
                       after = allAppInd[funi])
  }


  return(linesMat)
}
