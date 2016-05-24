merge.matrix <- function(
  x,
  y,
  ...
  ) {
  cn1 = colnames(x)
  cn2 = colnames(y)
  n1 = nrow(x)
  n2 = nrow(y)
  cnAll = unique(c(cn1, cn2))
  
  alldata = matrix(NA, n1 + n2, length(cnAll))
  colnames(alldata) = cnAll
  alldata[1:n1, match(cn1, cnAll)] = x
  alldata[n1 + (1:n2), match(cn2, cnAll)] = y
  alldata
}
