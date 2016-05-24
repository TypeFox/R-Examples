
# a JavaScript-like eval function
eval.js <- function(expr, envir=parent.frame(), enclos=if(is.list(envir)||is.pairlist(envir)) parent.frame())
{
  if (typeof(expr) != "character")
    return(expr)

  expr <- parse(text=expr)
  eval(expr, envir, enclos)
}

