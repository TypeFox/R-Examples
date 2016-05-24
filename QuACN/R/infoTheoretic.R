.infoTheoretic <- function(fvi, lambda=1000, nan.message="check your parameters") {
  fvi.sum <- sum(fvi)
  pis <- fvi/fvi.sum
  result <- list()
  result[["entropy"]] <- (-sum(pis*log2(pis)))
  if(is.nan(result[["entropy"]]))
    warning(paste("Entropy returned NaN (not a number):", nan.message))
  result[["distance"]] <- (lambda*(log2(length(pis)) - result[["entropy"]]))
  result[["pis"]] <- pis
  result[["fvis"]] <- fvi

  result
}

.infoTheoreticCoeff <- function(name, custom, diameter) {
  if (name == "lin")
    diameter:1
  else if (name == "quad")
    (diameter:1)^2
  else if (name == "exp")
    diameter * exp(0:(1 - diameter))
  else if (name == "const")
    rep(1, diameter)
  else if (name == "cust")
    custom
  else
    stop(paste("unknown coefficient setting '", name,
        "'; must be one of 'lin', 'quad', 'exp', 'const' or 'cust",
        sep=""))
}
