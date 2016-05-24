x2int <- function(x) 
{
  warn <- getOption("warn")
  options(warn = -1)
  rval <- as.integer(as.numeric(as.character(x)))
  options("warn" = warn)
  rval
}

