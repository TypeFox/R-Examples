se <-
function(x,
               na.rm = FALSE)
{
  if(!is.vector(x))
    stop("The function is only defined for vectors")
  if(na.rm)
    x <- x[!is.na(x)]
  sqrt(var(x) / length(x))
}
