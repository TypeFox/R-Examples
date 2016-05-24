imputeRough <- function(data, Init = "mean") {
  Pre <- my.dummy.df(data)
  if(Init == "mean") {
    output <- as.data.frame(apply(Pre, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)))
  } else {
    output <- as.data.frame(apply(Pre, 2, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x)))
  }
  Preds <- output[sapply(Pre, function(x) is.na(x))]
  list(Initials = Preds, Pre.Imputed = Pre, Imputed.Dataframe = output)
}