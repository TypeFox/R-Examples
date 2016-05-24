# Export objects with names in "toExport" from "envir" enviroment to "cluster".
# If ALL == TRUE it exports all the objects in "envir"
.clusterExport <- function(cluster, envir = parent.frame(), toExport = c(), ALL = FALSE)
{
  allNames <- toExport
  
  if(ALL) allNames <- c(allNames, ls(envir = envir))
  
  clusterExport(cluster, varlist = allNames, envir = envir)
}