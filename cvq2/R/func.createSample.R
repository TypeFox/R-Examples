func.createSample <-
function( input ){
  n <- nrow(input$modelData)
  #DEFAULT: Leave-One-Out cross validation, number of groups is equal to the number of elements 
  if( input$nFold == n )
    return( 1:n )
  else
    return( sample(1:n) )
}

