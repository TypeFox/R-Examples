func.output.regressionFormulaWithCoefficients <-
function( output, parameter ){
  if( is.null(output$writeTarget) )
    return()

  writeLines( "Regression formula with fitted parameters: ", con = output$writeTarget)

  start = 1
  param.i = 1
  text <- NULL

  if(names(output$coefficients)[1] == "(Intercept)")
    increment(start)

  for(i in start:NROW(output$coefficients) ){
    text[NROW(text)+1] <- paste(round(output$coefficients[[i]], output$round),"*",parameter[param.i])
    increment(param.i)
  }

  if( start == 2 )
    text[NROW(text)+1] <- round(output$coefficients[[1]], output$round)

  #use collapse instead sep, to combine the text vector
  cat(parameter[param.i],"=", paste(text,collapse=" + "), "\n", file = output$writeTarget)
}

