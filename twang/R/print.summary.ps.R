# Produces a summary table for ps object 
print.summary.ps <- function(x, ...)
{
      dots <- list(...)
      if(!is.null(dots$digits))
      obj <- round(x, digits = digits)
      else
      obj <- x
      class(obj) <- "matrix"
      print(obj)
      }

