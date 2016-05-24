process.histogram.variable <-
  function(variableName, concept, sym.obj.names) {
    suppressWarnings(sqldf(paste0("CREATE INDEX IF NOT EXISTS main.", variableName, " ON dataTable (", variableName, ")")))
    
    conceptColumns <- paste(concept, collapse = ", ")
    conceptConcatenation <- paste(concept, collapse = "||'.'||")
    categories <- sqldf(paste0("SELECT DISTINCT ", variableName, " FROM main.dataTable ORDER BY ", variableName))[[1]]
    
    result <- data.frame(rep("$H", length(sym.obj.names)), length(categories), check.names = F)
    colnames(result) <- c("$H", substr(variableName, 2, nchar(variableName) - 1))
    
     for (i in seq(from = 1, to = length(categories), by = 32)) {
      if (length(categories) - i + 1 >= 32)
        categoryGroup <- categories[i:(i+31)]
      else
        categoryGroup <- categories[i:length(categories)]
      
      queries <- character()
      for (category in categoryGroup) {
        queries <- c(queries, paste0("(SELECT SymObjNames, ifnull(freq, 0) AS '", category, "' FROM (SELECT SymObjNames FROM main.symObjTable) LEFT JOIN (SELECT ", conceptConcatenation, " AS concept, COUNT(", variableName, ") AS freq FROM main.dataTable WHERE ", variableName, " = '", category, "' GROUP BY ", conceptColumns, ") ON concept = SymObjNames)"))
      }
      queries <- paste(queries, collapse = " NATURAL JOIN ")
      result <- cbind(result, sqldf(paste0("SELECT * FROM ", queries))[-1])
    }
    
    totalFrequency <- sqldf(paste0("SELECT COUNT(", variableName, ") FROM main.dataTable GROUP BY ", conceptColumns, " ORDER BY ", conceptColumns))[[1]]
    
    for (j in 1:length(categories)) {
      result[j + 2] <- round(result[[j + 2]] / totalFrequency, 3)
    }
    
    return (result)
  }