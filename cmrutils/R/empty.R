`empty` <-
function(x)
{
  if(!is.null(x))
  {
    result <- x
    result[] <- NA
  }
  else
    result <- NULL
  
  return(result)
}
