parse.random.bayesx <-
function(term, data)
{
  ra <- eval(parse(text = term))
  term <- ra$term
  by <- ra$by
  h.random <- NULL
  if(!is.null(ra$ins)) {
    if(is.null(ra$data))
      ra$data <- data
    h.random <- parse.bayesx.input(ra$formula, ra$data, 
      ra$weights, ra$subset, ra$offset, 
      ra$na.action, ra$contrasts, ra$control)
    h.random$family <- "gaussian_re"
  }    

  return(list(h.random = h.random, term = term, by = by))
}

