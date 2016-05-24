safe_efficiency <- function(fit, type) {
  res <- try(efficiency(fit, type = type, plot = FALSE)[c(type, "eff", "fluo")], 
             silent = TRUE)
  if (class(res) == "try-error") {
    res <- rep(NaN, length(c(type, "eff", "fluo"))) 
  } else {
    if (length(res[["eff"]]) > 1)
      res$eff <- NaN
  }
  unlist(res)
}