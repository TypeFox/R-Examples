help.map.name <- function(x)
{
  if(is.null(x))
    return("")
  x <- splitme(x)
  if(resplit(x[1L:2L]) == "s(") {
    x <- resplit(c("sfoofun", x[2L:length(x)]))
    x <- eval(parse(text = x), envir = parent.frame())
    if(is.null(x))
      x <- "bayesxmap"
  } else  x <- "bayesxmap"

  return(x)
}

