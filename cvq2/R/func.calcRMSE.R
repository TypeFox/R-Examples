func.calcRMSE <-
function( yiCalc, yi, nu=0 ){
  n <- NROW(yiCalc) - nu
  
  return( sqrt( sum( (yiCalc - yi) * (yiCalc - yi) ) / n ) )
}

