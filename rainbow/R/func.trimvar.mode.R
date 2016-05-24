`func.trimvar.mode` <- function(data, trim = 0.25)
{
  functions = t(data$y)
  lista = depth.mode(data, trim = trim)$ltrim; 
  func.var(functions[lista, ])
}

