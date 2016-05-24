summary.bayesx <-
function(object, model = NULL, 
  digits = max(3, getOption("digits") - 3), ...)
{
  res <- get.model(object, model)
  attr(res, "digits") <- digits
  attr(res, "args") <- list(...)
  class(res) <- "summary.bayesx"

  return(res)
}

