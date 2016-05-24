getIndex <-
function(idx){
  rrow <- ceiling((-1 + sqrt(8*idx + 1))/2)
  ccol <- idx - rrow*(rrow - 1)/2
  ans <- c(rrow, ccol)
  return(ans)
}
