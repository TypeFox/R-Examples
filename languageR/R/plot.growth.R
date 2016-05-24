`plot.growth` <-
function(x, w = "all", ...) {
  if (!is(x, "growth")) stop("argument should be a growth object")
  res = x@data$data
  if (w == "all") {
    par(mfrow = c(2, 4))
    plot(res$Tokens, res$Types, xlab = "tokens", ylab = "types")
    plot(res$Tokens, res$HapaxLegomena/res$Tokens, xlab = "tokens",
      ylab = "growth rate")
    plot(res$Tokens, res$TypeTokenRatio, xlab = "tokens", ylab = 
      "type-token ratio")
    plot(res$Tokens, res$Lognormal, xlab = "tokens", 
      ylab = "mean log frequency")
    plot(res$Tokens, res$Herdan, xlab = "tokens", ylab = "Herdan's C")
    plot(res$Tokens, res$Guiraud, xlab = "tokens", ylab = "Guiraud's R")
    plot(res$Tokens, res$Yule, xlab = "tokens", ylab = "Yule's K")
    plot(res$Tokens, res$Zipf, xlab = "tokens", ylab = "Zipf slope")
    par(mfrow = c(1, 1))
  } else {
    if (w %in% colnames(res)) {
      plot(res$Tokens, res[, w], xlab="Tokens", ylab=w)
    } else
      stop("w should specify a valid column name")
  }
}

