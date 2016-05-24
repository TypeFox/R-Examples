countnotna <-
function(x)
# countnotna()  - find no of non-na elements in a vector
{
  l <- length(x[!is.na(x)])
  return(l)
}
