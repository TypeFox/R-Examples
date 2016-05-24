func.output.coefficients <-
function( output ){
  if( is.null(output$writeTarget) )
    return()

  func.output.linearFormula( output )

  intercept <- 0

  # output coefficents, formatted
  # NROW gross, da Vektor
  for(i in 1:NROW(output$coefficients) ){
    param <- letters[i - intercept]
    value <- round(output$coefficients[i], output$round)
    
    if(names(output$coefficients)[i] == "(Intercept)"){
      param <- "const"
      increment(intercept)
    }
  
    writeLines( paste( param,":", value), con = output$writeTarget )
  }
}

