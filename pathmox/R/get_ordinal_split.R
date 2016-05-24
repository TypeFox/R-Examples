#' @title Binary partitions of an ordinal variable
#' 
#' @description
#' Internal function. \code{get_ordinal_split} is called by \code{pathmox}
#' 
#' @param cat.exe vector of categories
#' @export
#' @keywords internal
get_ordinal_split <- function(cat.exe)
{
  a = length(cat.exe)
  if (a == 2)  # binary variable
  {
    part = cat.exe[1]
    antipar = cat.exe[2]
  }
  if (a > 2)   # more than 2 categories
  {
    part = as.list(1:(a-1))
    antipar = part
    for (i in 1:(a-1))
    {
      part[[i]] = cat.exe[1:i]
      antipar[[i]] = setdiff(cat.exe, part[[i]])
    }
  }
  # output
  list(par1 = part, par2 = antipar)
}
