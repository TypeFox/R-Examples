`transforming.fnc` <-
function(y, fun) {
  if (is.function(fun)) {
    return(fun(y))
    #if (fun == "exp") {
    #  return(exp(y))
    #} else {
    #  return(plogis(y))
    #}
  } else 
    return(y)
}

