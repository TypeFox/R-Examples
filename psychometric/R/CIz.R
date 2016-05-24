"CIz" <-
function (z, n, level=.95) 
 { 
  noma <- 1-level
  sez <- SEz(n)
 zs <- - qnorm(noma/2)
 mez <- zs*sez
 lcl <- z - mez
 ucl <- z + mez
 mat <- list(lcl, ucl)
 return(as.numeric(mat))
 }

