desc <- function(table, tolower=TRUE, dots=FALSE, ...)
{
  on.exit(suppressWarnings(dbUnloadDriver(dbDriver("Oracle"))))

  ## 1  Fetch table description
  query <- paste("SELECT * FROM", table)
  output <- dbColumnInfo(dbSendQuery(dbConnect(dbDriver("Oracle"),...), query))

  ## 2  Handle row and column names
  attr(output,"row.names") <- seq_len(nrow(output))
  if(tolower)
  {
    output$name <- tolower(output$name)
    output$type <- tolower(output$type)
  }
  if(dots)
    output$name <- chartr("_", ".", output$name)

  ## 3  Add row count info
  splitname <- toupper(unlist(strsplit(table, "\\.")))
  select.from <- "SELECT num_rows,last_analyzed FROM all_tables"
  if(length(splitname) == 1)
    where <- paste("WHERE table_name='", splitname, "'", sep="")
  else
    where <- paste("WHERE owner='", splitname[1], "' AND table_name='", splitname[2], "'", sep="")
  query <- paste(select.from, where)
  rows.date <- dbGetQuery(dbConnect(dbDriver("Oracle"),...), query)
  attr(output, "rows") <- if(nrow(rows.date)==1) rows.date$NUM_ROWS else as.numeric(NA)
  attr(output, "analyzed") <- if(nrow(rows.date)==1) rows.date$LAST_ANALYZED else as.character(NA)

  ## 4  Show on screen
  cat("\n")
  print(output, right=FALSE)
  cat("\n")
  if(nrow(rows.date) > 0)
  {
    indent <- rep(" ", nchar(nrow(output))+1)
    cat(indent, attr(output,"rows"), " rows on ", as.character(attr(output,"analyzed")), "\n\n", sep="")
  }

  invisible(output)
}
