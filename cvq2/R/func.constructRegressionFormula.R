func.constructRegressionFormula <-
function( col.name ){
  #y column
  col.formula = paste( col.name[NROW(col.name)], "~" )

  for( i in 1:(NROW(col.name)-1) ){
    tmp.sign = "+"
    if( i == 1 )
      tmp.sign = ""

    col.formula = paste( col.formula, tmp.sign, col.name[i] )
  }
  return( as.formula(col.formula) )
}

