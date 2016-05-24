get.s1 <-
function( i, block.obj, Xa = NA, Y = NA, r2.cut, s1high, s1low, max.iter, tol ){
  s1 <- NA
  while( is.na(s1) ){
    if(isS4(block.obj)){
      s1 <- p0tos1( p0 = i, block.obj, r2.cut = r2.cut, s1high = s1high, s1low=s1low, max.iter = max.iter, tol = tol )
    }else{
      s1 <- p0tos1( p0 = i, block.obj = NA, Xa = Xa, Y = Y, r2.cut = r2.cut, s1high = s1high, s1low=s1low, max.iter = max.iter, tol = tol )
    }
    i <- i+1
  }
  return(s1)
}

