GERule.Concat <- function(x,y) {
  lx = length(x)
  ly = length(y)
  for (i in 1:ly) {
    x[lx + i] = y[i]
  }
  
  x
}

`+.GERule` <- GERule.Concat
c.GERule <- function(..., recursive = FALSE) {
  l = list(...)
  
  ret = l[[1]]
  if (length(l) > 1) {
    for (i in 2:length(l)) {
      ret = GERule.Concat(ret, l[[i]])
    }
  }
  
  ret
}

