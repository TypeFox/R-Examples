readDescriptionSDML <- function(x)
{
  if(is.null(x)) return(NULL)

  list(
       title      = xmlValue(x[["title"]]),
       source     = xmlValue(x[["source"]]),
       date       = xmlValue(x[["date"]]),
       version    = xmlValue(x[["version"]]),
       comment    = xmlValue(x[["comment"]]),
       creator    = xmlValue(x[["creator"]]),
       properties = readProperties(x[["properties"]])
       )
}
