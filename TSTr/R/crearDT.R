crearDT <-
  function(col1 = character(), 
           col2 = character()) {
    DT <- data.table(V1 = col1, V2 = col2)
    return(DT)
  }
