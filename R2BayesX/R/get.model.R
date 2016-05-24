get.model <- function(x, model)
{
  elmts <- c("formula", "bayesx.setup", "bayesx.prg", "bayesx.run",   
    "fitted.values", "residuals", "response", "effects", "smooth.hyp", 
    "random.hyp", "fixed.effects", "variance", "model.fit", "call")
  if(!any(names(x) %in% elmts)) {
    if(!is.null(model)) {
      if(is.character(model)) {
        if(all(is.na(model <- pmatch(model, names(x)))))
          stop("argument model is specified wrong!")
      } else {
        if(max(model) > length(x) || is.na(model) || min(model) < 1) 
          stop("argument model is specified wrong!")
      }
      x <- x[model]
    }
  } else x <- list(x)
  class(x) <- "bayesx"

  return(x)
}

