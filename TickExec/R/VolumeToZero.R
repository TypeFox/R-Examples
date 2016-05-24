

VolumeToZero <- function (df) {
  df[which(is.na(df) == TRUE, arr.ind = TRUE)] = 0
  return (df)
}