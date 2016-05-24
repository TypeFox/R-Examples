cprob <-
function(object, model = NULL, term = NULL, ...)
{
  if(is.null(term))
    term <- 1L
  object <- get.model(object, model)
  if(length(object) > 1L) {
    rval <- list()
    for(i in 1L:length(object))
      rval[[i]] <- attr(object[[i]]$effects[[term]], "contourprob")
    } else rval <- attr(object[[1L]]$effects[[term]], "contourprob")
  if(is.null(rval) || length(rval) < 1L)
    rval <- NA

  return(rval)
}

