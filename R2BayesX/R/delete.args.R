delete.args <- function(fun = NULL, args = NULL, not = NULL, package = NULL)
{
  if(is.character(fun) & !is.null(package))
    fun <- eval(parse(text = paste(package, paste(rep(":", 3), collapse = ""), fun, sep = "")))
  nf <- names(formals(fun))
  na <- names(args)
  for(elmt in na)
    if(!elmt %in% nf) {
      if(!is.null(not)) {
        if(!elmt %in% not)
          args[elmt] <- NULL
      } else args[elmt] <- NULL
    }

  return(args)
}

