delete.NULLs <- function(x.list) 
{
  x.list[unlist(lapply(x.list, length) != 0)]
}
