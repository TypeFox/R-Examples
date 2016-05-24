func.sortDataColumns <-
function( input, sortData, writeOutputTarget){

  new.data <- cbind(sortData[,all.vars(input$regressionFormula)[2]])
  colnames(new.data) <- all.vars(input$regressionFormula)[2]

  if( NCOL(sortData) > 2 ){
    for( i in 3:NCOL(sortData)){
      new.data <- cbind(new.data, sortData[,all.vars(input$regressionFormula)[i]] )
      colnames(new.data)[i-1] <- all.vars(input$regressionFormula)[i]
    }
  }

  # y as last column
  new.data <- cbind(new.data, sortData[,all.vars(input$regressionFormula)[1]])
  colnames(new.data)[NCOL(sortData)] <- all.vars(input$regressionFormula)[1]

  rownames(new.data) <- 1:NROW(new.data)
#  writeLines (colnames(input$modelData) == colnames(new.data), con = writeOutputTarget )

  tmp.diff <- TRUE
  for( i in 1:NCOL(sortData)){
    if( colnames(sortData)[i] != colnames(new.data)[i] ){
      tmp.diff <- FALSE
    }
  }
  if( !tmp.diff && !is.null(writeOutputTarget) ){
    writeLines("WARNING, the columns of your input table have to be reordered to fit the given regression formula", con = writeOutputTarget )
    writeLines("Input of data:", con = writeOutputTarget )
    write.table(sortData, file = writeOutputTarget, sep="\t", row.names = FALSE )
    writeLines( paste("New Sort of columns according to formula:",deparse(input$regressionFormula)), con = writeOutputTarget )
    write.table(new.data, file = writeOutputTarget, sep="\t", row.names = FALSE )
    writeLines("", con = writeOutputTarget )
  }
  return( as.data.frame(new.data) )
}

