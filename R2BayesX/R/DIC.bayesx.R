DIC.bayesx <- function(object, ...)
{
  object <- c(object, ...)
  val <- extract.model.diagnostic(object, 1L:length(object), "DIC", FALSE)
  val <- data.frame(pd = extract.model.diagnostic(object, 1L:length(object), "pd", FALSE), 
    DIC = val)
  Call <- match.call()
  Call$k <- NULL
  row.names(val) <- if(nrow(val) > 1) as.character(Call[-1L]) else ""

  return(val)
}

