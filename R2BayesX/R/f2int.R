f2int <- function(x, type = 1L) 
{
  if(!is.factor(x))
    x <- as.factor(x)
  if(type != 3L) {
    levels(x) <- nl <- 1:nlevels(x)
    x <- factor(x, levels = nl, labels = nl)
    x <- as.integer(x)
    if(type != 2L)
      x <- x - 1L
    if(min(x, na.rm = TRUE) < 0)
      x <- x + abs(min(x, na.rm = TRUE))
  } else {
    warn <- getOption("warn")
    options("warn" = -1)
    x2 <- as.integer(as.character(x))
    x <- if(all(is.na(x2))) as.integer(x) else x2
    options("warn" = warn)
  }

  return(x)
}
