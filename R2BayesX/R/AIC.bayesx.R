AIC.bayesx <- function(object, ..., k = NULL)
{
  object <- c(object, ...)
  val <- extract.model.diagnostic(object, 1L:length(object), "AIC", FALSE)
  val <- data.frame(df = extract.model.diagnostic(object, 1L:length(object), "df", FALSE), 
    AIC = val)

  Call <- match.call()
  Call$k <- NULL
  row.names(val) <- if(nrow(val) > 1) as.character(Call[-1L]) else ""

  return(val)
}

