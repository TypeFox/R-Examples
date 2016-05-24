gls.approx.logistic <- function(snpdata, leftvar, rightvars = NULL, outvar = paste(leftvar, "star", sep = ""), weightvar = "weight") {
  stopifnot(is.data.frame(snpdata$data))
  stopifnot(leftvar %in% names(snpdata$data))
  stopifnot(all(rightvars %in% names(snpdata$data)))
  stopifnot(all(snpdata$data[[leftvar]] %in% c(0, 1, NA)))
  if (outvar %in% names(snpdata$data)) warning("overwriting ", outvar, " in snpdata$data")
  if (weightvar %in% names(snpdata$data)) warning("overwriting ", weightvar, " in snpdata$data")
  null <- as.formula(paste(leftvar, paste(c("1", rightvars), collapse = "+"), sep = "~"))
  print(null)
  p <- 1/(1 + exp(-predict(glm(null, family = "binomial", data = snpdata$data, na.action = na.exclude))))
  snpdata$data[[weightvar]] <- p*(1-p)
  snpdata$data[[outvar]] <- (snpdata$data[[leftvar]] - p)/snpdata$data[[weightvar]]
  return(snpdata)
}
