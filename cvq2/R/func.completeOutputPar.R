func.completeOutputPar <-
function( output, extOut, extOutFile ){
  #change DEFAULT
  if( !is.null(extOutFile) )
    extOut <- TRUE

  if( extOut ){
    if( is.null(extOutFile) )
      output$writeTarget = stdout()
    #redirect output to file
    else{
      output$toFile <- TRUE
      output$writeTarget = file(extOutFile, open = "w")
    }
  }

  return( output )
}

