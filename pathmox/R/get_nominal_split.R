#' @title Binary partitions of a nominal variable
#' 
#' @description
#' Internal function. \code{get_nominal_split} is called by \code{pathmox}
#' 
#' @param v vector containing the nominal categories
#' @export
#' @keywords internal
get_nominal_split <- function(v)
{
  # length of vector
  k = length(v)
  # number of binary partitions
  num.parts = (2^(k-1)) - 1   
  # stop limit of order of partitions
  if ((k%%2) == 0) stop_parts = k/2 else stop_parts = (k-1)/2
  # list for storing the partitions
  part = as.list(1:num.parts)
  # complementary list for part
  antipar = part   
  
  # FUNCTION(S)
  # internal function for evaluating candidate splits
  eval.cand = function(x, y, aux)
  {
    # x: list of partitions
    # y: candidate split to be included in part
    # aux: num of parts already included in part
    acum = 0
    for (i in 1:(aux-1))
    {
      true.false = setequal(x[[i]], y)
      if (true.false) acum = acum + 1 else acum = acum + 0
    }
    return(acum)
  }
  
  # first order partitions
  aux = 1
  for (i in 1:k)
  {
    part[[aux]] = v[i]
    antipar[[aux]] = setdiff(v, part[[i]])
    aux = aux + 1
  }
  
  # partitions of order greater than one
  for (ord in 2:stop_parts)
  {
    if (stop_parts < ord) break
    # partitions of the immediate past order
    p1 = aux - choose(k, (ord-1))
    p2 = aux - 1
    part.ant = part[p1:p2]
    antipar.ant = antipar[p1:p2]    
    for (i in 1:length(part.ant))
    { 
      part1 = part.ant[[i]]
      anti1 = antipar.ant[[i]]
      for (j in 1:length(anti1))
      {
        candidat = c(part1, anti1[j])
        cand = eval.cand(part, candidat, aux)
        if (cand == 0)
        {
          part[[aux]] = candidat
          antipar[[aux]] = setdiff(v, candidat)
          aux = aux + 1
        } else
        {
          next
        }
      }
    }
  }    
  # in case that k is even, only select partitions until num.parts
  if (length(part) > num.parts)
  {
    part = part[1:num.parts]
    antipar = antipar[1:num.parts]
  }   
  # output
  list(par1 = part, par2 = antipar)
}
