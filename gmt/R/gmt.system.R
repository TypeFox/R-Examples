gmt.system <- function(cmd, file=NULL, append=FALSE)
{
  output <- system(cmd, intern=TRUE)
  open <- if(append) "a" else "w"  # append or write

  if(is.character(file))
  {
    con <- file(file, open)
    writeLines(output, con)
    close.connection(con)
  }

  invisible(output)
}
