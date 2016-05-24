term.names <- function(x)
{
  n <- length(x)
  names <- vector("list", n)
  for(i in 1L:length(x)) {
    if(is.sm(x[i]))
      names[[i]] <- eval(parse(text = x[i]))$term
    else
      names[[i]] <- x[i]
  }
  
  return(names)
}

