r2gmt <- function(x, datafile, append=FALSE)
{
  is.digit <- function(x) suppressWarnings(!is.na(as.numeric(x)))

  ## 1  Import x if it's a file
  if(is.character(x) && file.exists(x))
  {
    char <- substring(scan(x,what="",n=1,quiet=TRUE), 1, 1)  # assume 1 line header if first non-whitespace char ...
    header <- all(char!="#", char!="-", char!=".", !is.digit(char))  # ... is not comment, minus, point, or number
    sep <- if(length(grep(",",readLines(x,n=1)))>1) "," else ""  # assume comma separator if first line has comma
    x <- read.table(x, sep=sep, header=header)
  }

  ## 2  Now that x is a data object, export it
  if(append)
    write(">", datafile, append=TRUE)
  write.table(x, datafile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=append)

  invisible(x)
}
