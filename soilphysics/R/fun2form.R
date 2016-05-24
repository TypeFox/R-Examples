fun2form <-
function(fun, y = NULL)
{
   if (!inherits(fun, "function"))
       stop("'fun' must be an object of class 'function'!")
   if (!is.null(y) & !inherits(y, "character"))
       stop("'y' must be a 'character' which is going to define the left side of the formula!")
   form <- as.formula(paste(y, "~", deparse(fun)[2]))
   return(form)
}
