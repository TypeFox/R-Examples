tables <- function(owner="%", table="%", space="%", tolower=TRUE, ...)
{
  ## 1  Prepare query
  select.from <- "SELECT owner,table_name,tablespace_name,num_rows,last_analyzed FROM all_tables"
  where <- paste("WHERE owner LIKE '", toupper(owner), "' ESCAPE '\\' ",
                 "AND table_name LIKE '", toupper(table), "' ESCAPE '\\' ",
                 "AND tablespace_name LIKE '", toupper(space), "' ESCAPE '\\'", sep="")
  query <- paste(select.from, where)

  ## 2  Run query
  output <- sql(query, ...)

  ## 3  Format output
  names(output) <- c("owner", "table", "space", "rows", "analyzed")
  if(tolower)
    output[1:3] <- sapply(output[1:3], tolower)

  ## 4  Show on screen
  cat("\n")
  print(output, right=FALSE)
  cat("\n")

  invisible(output)
}
