func.calcXSquare <-
function( yiCalc, yi, yMean ){
  return( 
    1 -
    sum( (yiCalc - yi) * (yiCalc - yi) ) /
    sum( (yi - yMean) * (yi - yMean) )
  )
}

