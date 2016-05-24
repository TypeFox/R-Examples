classic.to.sym <-
  function(dataTable, concept, variables, variables.types) {
    
    if (length(variables) != length(variables.types)) {
      stop("variables and variables.types must have the same length")
    }
    
    dataTable <- dataTable[, which(colnames(dataTable) %in% c(variables, concept))]
    
    concept <- paste0("[", concept, "]")
    variables <- paste0("[", variables, "]")
    conceptColumns <- paste(concept, collapse = ", ")
    
    sqldf() 
    
    sqldf(paste0("CREATE INDEX main.concept_index ON dataTable (", conceptColumns, ")"))
    
    sym.obj <- sqldf(paste0("SELECT DISTINCT ", conceptColumns, " FROM main.dataTable ORDER BY ", conceptColumns))
    sym.obj.names <- do.call("paste", args = c(sym.obj, sep = "."))
    
    symObjTable <- data.frame(SymObjNames = sym.obj.names)
    
    sqldf("SELECT SymObjNames FROM symObjTable")
    
    meta.data <- list()
    for (i in 1:length(variables)) {
      
      switch (variables.types[[i]],
              '$C' = {
                meta.data[[i]] <- process.continuum.variable(variables[[i]], conceptColumns)
              },
              '$I' = {
                meta.data[[i]] <- process.interval.variable(variables[[i]], conceptColumns)
              },
              '$H' = {
                meta.data[[i]] <- process.histogram.variable(variables[[i]], concept, sym.obj.names)
              },
              '$S' = {
                meta.data[[i]] <- process.set.variable(variables[[i]], concept, sym.obj.names)
              },
              stop("Invalid variable type"))
    }
    
    sqldf()  
    
    meta.data <- data.frame(meta.data, check.names = F)
    rownames(meta.data) <- sym.obj.names
    
    colnames(meta.data)[colnames(meta.data) == "'$C'"] <- "$C"
    colnames(meta.data)[colnames(meta.data) == "'$I'"] <- "$I"
    
    return (newSobject(meta.data))
  }