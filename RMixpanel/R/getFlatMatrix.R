getFlatMatrix <- function(
  df
  ) {
### MP, 2015
  res = matrix("", nrow(df), ncol(df))
  for(i in 1:ncol(df)) {
    dfi = unlist(lapply(df[, i], paste, collapse=","))
    res[, i] = dfi
  }
  colnames(res) = colnames(df)
  res
}
