process.continuum.variable <-
function(variableName, conceptColumns) {
  varType <- "$C"
  return (sqldf(paste0("SELECT '", varType, "', round(AVG(", variableName, "), 2) AS ", variableName, " FROM main.dataTable GROUP BY ", conceptColumns, " ORDER BY ", conceptColumns)))
}
