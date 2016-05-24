big.peaks.col <-
function(x, tre){
  r <- rle(x)
  v <- which(rep(x = diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2, times = r$lengths))
  pos <- v[which(x[v] > tre)] #position of the real peaks
  hei <- x[pos]
  out <- list(pos=pos,hei=hei)
  return(out)
}
