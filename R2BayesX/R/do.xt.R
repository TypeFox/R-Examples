do.xt <- function(term, object, not = NULL, noco = FALSE)
{
  if(!is.null(object$xt)) {
    names.xt <- names(object$xt)
    if(is.null(not))
      not <- "not"
    count <- 1
    co <- ","
    if(noco)
      co <- NULL
    for(name in names.xt) {
      if(count > 1)
        co <- ","
      if(!name %in% not) {
        if(name %in% c("full", "catspecific", "center", "derivative", "nofixed") || 
          is.logical(object$xt[[name]])) {
          if(is.logical(object$xt[[name]])) {
            if(object$xt[[name]])
              term <- paste(term, co, name, sep = "")
          } else term <- paste(term, co, name, sep = "")
          count <- count + 1
        } else {
          term <- paste(term, co, name, "=", object$xt[name], sep = "")
          count <- count + 1
        }
      }
    }
  }

  return(term)
}

