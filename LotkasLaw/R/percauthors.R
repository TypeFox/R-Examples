percauthors <- function(Table, Sums)
{
  justy <- Table[,2:2]
  newcol <- justy/Sums[2]
  return(newcol)
  
}
