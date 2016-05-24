GCV.bayesx <- function(object, ...)
{
  object <- c(object, ...)
  val <- extract.model.diagnostic(object, 1L:length(object), "GCV", FALSE)
  val <- data.frame(df = extract.model.diagnostic(object, 1L:length(object), "df", FALSE), 
    GCV = val)

  Call <- match.call()
  Call$k <- NULL
  row.names(val) <- if(nrow(val) > 1) as.character(Call[-1L]) else ""

  return(val)
}

