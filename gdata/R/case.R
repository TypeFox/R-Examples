case <- function(x, ..., default=NA)
{
  magic <- "....default...."
  alternatives <- c(...,"....default...."=magic)

  x <- as.character(x)
  retval <- factor(
                   x,
                   levels=alternatives,
                   labels=names(alternatives)
                   )
  levels(retval)[length(alternatives)] <- as.character(default)  
  retval[is.na(retval) & !is.na(x)] <- default

  retval
}
