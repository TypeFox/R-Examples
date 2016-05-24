`getRaggedSeq` <-
function(data) {
  if( dim(data)[2]!=2)
    warning("data should have 2 columns")
  freq = tapply(data[,2], as.character(data[,1]), function(x) length(unique(x)))  
  result = cumsum(freq)
  result = c(start=1, 1+result)
  names(result) = c(names(result)[-1], "end")
  # get rid of trailing white space
  names(result) = gsub("[[:space:]]+$", "", names(result))
  return(result)
}

