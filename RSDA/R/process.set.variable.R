process.set.variable <-
  function(variableName, concept, sym.obj.names) {
    suppressWarnings(sqldf(paste0("CREATE INDEX IF NOT EXISTS main.", variableName, " ON dataTable (", variableName, ")")))
    
    conceptColumns <- paste(concept, collapse = ", ")
    conceptConcatenation <- paste(concept, collapse = "||'.'||")
    categories <- sqldf(paste0("SELECT DISTINCT ", variableName, " FROM main.dataTable ORDER BY ", variableName))[[1]]
    
    result <- data.frame(rep("$S", length(sym.obj.names)), length(categories), check.names = F)
    colnames(result) <- c("$S", substr(variableName, 2, nchar(variableName) - 1))
    
     for (i in seq(from = 1, to = length(categories), by = 64)) {
      if (length(categories) - i + 1 >= 64)
        categoryGroup <- categories[i:(i+63)]
      else
        categoryGroup <- categories[i:length(categories)]
      
      queries <- character()
      for (category in categoryGroup) {
        queries <- c(queries, paste0("(SELECT SymObjNames, SymObjNames IN (SELECT DISTINCT ", conceptConcatenation, " FROM main.dataTable WHERE ", variableName, " = '", category, "') AS '", category, "' FROM main.symObjTable)"))
      }
      queries <- paste(queries, collapse = " NATURAL JOIN ")
      result <- cbind(result, sqldf(paste0("SELECT * FROM ", queries))[-1])
    }
    
    return (result)
  }
