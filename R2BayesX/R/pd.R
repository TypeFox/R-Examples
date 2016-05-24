pd <-
function(object, model = NULL, print.names = FALSE)
{
  return(extract.model.diagnostic(object, model, "pd", print.names))
}

