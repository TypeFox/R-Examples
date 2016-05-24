

PriceToNA <- function (df) {
  df[which(df == 0, arr.ind = TRUE)] = NA
  return (df)
}