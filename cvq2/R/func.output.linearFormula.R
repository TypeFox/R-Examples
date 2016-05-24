func.output.linearFormula <-
function( output ){
  if( is.null(output$writeTarget) )
    return()

  coeffLast <- "z*Xn"
  if(names(output$coefficients)[1] == "(Intercept)")
    coeffLast <- "const"

  writeLines( paste("General Equation: Y = a*X1 + b*X2 + ... + ", coeffLast) , con = output$writeTarget )
}

