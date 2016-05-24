views <- function(owner="%", view="%", tolower=TRUE, ...)
{
  ## 1  Prepare query
  select.from <- "SELECT owner,view_name FROM all_views"
  where <- paste("WHERE owner LIKE '", toupper(owner), "' ESCAPE '\\' ",
                 "AND view_name LIKE '", toupper(view), "' ESCAPE '\\'", sep="")
  query <- paste(select.from, where)

  ## 2  Run query
  output <- sql(query, ...)

  ## 3  Format output
  names(output) <- c("owner", "view")
  if(tolower)
    output[1:2] <- sapply(output[1:2], tolower)

  ## 4  Show on screen
  cat("\n")
  print(output, right=FALSE)
  cat("\n")

  invisible(output)
}
