df2m2 <-
function(x)
{
  x$intnr <- NULL
  x$paramnr <- NULL
  vars <- as.character(x$varname)
  x$varname <- NULL
  x <- as.matrix(x)
  vars[vars == "const"] <- "(Intercept)"
  rownames(x) <- vars

  return(x)
}

