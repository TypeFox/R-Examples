generatelens <- function(maxLen, minLen, stdevLen, nString, distLen){
  if(distLen == 'normal'){
    medLen <- (maxLen+minLen)/2
    len <- rtruncnorm(n = nString, a = minLen, b = maxLen, mean = medLen, sd = stdevLen)
    len <- round(x = len)
  }
  else if(distLen == 'uniform'){
    len <- sample(minLen:maxLen, nString, replace = TRUE)
  }
  else{
    len <- "err in generatelens"
  }
  len
}
