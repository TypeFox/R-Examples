pasteRound <- function (..., digits=16, sep=" ", collapse=NULL) {
  args <- list(...)
  if (length(args) == 0)
    if (length(collapse) == 0)
      character(0)
    else ""
  else{
    for(i in seq_along(args))
      if(is.numeric(args[[i]])) args[[i]] <-  as.character(round(args[[i]], digits))
      else args[[i]] <- as.character(args[[i]])
   paste(args, sep, collapse)
  }
}



