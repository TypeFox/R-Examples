mean.rivernet <- function(x,y=NA,na.rm=FALSE,...)
{
  rivernet <- x
  lengths  <- rivernet$attrib.reach$length
  orders   <- rivernet$attrib.reach$streamorder
  if ( length(y)==1 & is.na(y[1]) ) return(mean(lengths,na.rm=na.rm,...))
  if ( length(y) != length(lengths) ) return(NA)
  if ( na.rm )
  {
    ind <- !is.na(y) & !is.na(lengths) & !is.na(orders)
    if ( sum(ind) == 0 ) return(NA)
    return(sum(y[ind]*lengths[ind]*orders[ind])/sum(lengths[ind]*orders[ind]))
  }
  else
  {
    return(sum(y*lengths*orders)/sum(lengths*orders))
  }
}
