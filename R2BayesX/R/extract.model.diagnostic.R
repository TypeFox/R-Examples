extract.model.diagnostic <- function(object, model, what, print.names = FALSE)
{
  object <- get.model(object, model)
  dg <- df <- NULL
  for(i in 1L:length(object)) {
    tmp <- eval(parse(text = paste("object[[i]]$model.fit$", what, sep = "")))
    if(is.null(tmp))
      tmp <- NA
    if(what == "logLik")
      df <- c(df, object[[i]]$model.fit$df)
    dg <- c(dg, tmp)
  }
  if(!is.null(dg) && print.names)
    names(dg) <- names(object)
  if(any(is.na(dg)))
    warning(paste("information on", what, "is not available!"), call. = FALSE)
  if(what == "logLik") {
    names(df) <- names(object)
    attr(dg, "df") <- df
  }

  return(dg)
}

