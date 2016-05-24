`nan2na` <-
function(x)
{
  if(!is.null(x))
  {
    result <- x
    result[!is.finite(result)] <- NA
  }
  else
    result <- NULL
  
  return(result)
}
