vfsort <- function( vf, decreasing = FALSE ) {
  
# sort data by id, eye, date, and time
  vf <- vf[order( vf$ttime, decreasing = decreasing ),]
  vf <- vf[order( vf$tdate, decreasing = decreasing ),]
  vf <- vf[order( vf$seye,  decreasing = decreasing ),]
  vf <- vf[order( vf$id,    decreasing = decreasing ),]
  
  return( vf )
  
}