read.bayesx.output <-
function(dir, model.name = NULL)
{
  if(is.null(dir))
    stop("no directory specified!")
  if(is.null(model.name))
    bm <- search.bayesx.models(dir)
  else
    bm <- model.name
  k <- length(bm)
  rval <- list()
  for(j in 1:k) {
    cmd <- paste("rval$", bm[j], "<-read.bayesx.model.output('", dir, "','", bm[j], "')", sep = "")
    eval(parse(text = cmd))
  }
  if(length(rval) < 2L)
    rval <- rval[[1L]]
  class(rval) <- "bayesx"

  return(rval)
}

