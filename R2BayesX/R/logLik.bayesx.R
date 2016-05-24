logLik.bayesx <-
function(object, ...)
{
  args <- list(...)
  model <- args$model
  if(!is.null(args$print.names))
    print.names <- args$print.names
  else
    print.names <- FALSE 
  return(extract.model.diagnostic(object, model, "logLik", print.names))
}

