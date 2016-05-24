`func.trimvar.RPD` <- function(functions, trim = 0.25, deriv = c(0, 1))
{
  lista = depth.RPD(functions, trim = trim, deriv = deriv)$ltrim; 
  func.var(functions[lista,])
}
