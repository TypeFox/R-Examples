parConstrainFun <- function(x, parameters, full, fm2, comp, G)
  {
   vapply(parameters[min(x):max(x)], FUN = constrainFun, FUN.VALUE = vector("numeric", length = 1), full = full, fm2 = fm2, comp = comp, G = G)
 }

