#============================================================================================
# Compute STI values
#============================================================================================
sti <- function(values, scale) {
  if (scale <= 0) {
    stop('Error: scale must be a positive value bigger or equal to one.')
  }
  if (length(values) == 0) {
    stop('Error: check your input data, its length is zero. It must be a vector of numeric values.')
  }
  if (length(values) < 360) {
    cat('Warning: just to remind you than a time series longer than 30 years is recommended.\n')
  }
  values <- as.numeric(values)
  nb.data <- length(values)
  sti.values <- rep(NA, nb.data)
  ts <- zoo(values)
  weighted.ts <- c(rep(NA, times=scale-1))
  weighted.ts <- c(weighted.ts, as.vector(rollapply(ts, width=scale, align="right", by=1, FUN=mean, na.rm = F)))
  for (month in 1:12) {
    indices <- c(1:nb.data)[seq(month, length(c(1:nb.data)), 12)]
    valid.weighted.ts <- weighted.ts[seq(month, length(weighted.ts), 12)][!is.na(weighted.ts[seq(month, length(weighted.ts), 12)])]
    n <- fitdist(valid.weighted.ts, distr="norm", method="mle")
    sti.values[indices] <- (weighted.ts[indices] - n$estimate["mean"]) / n$estimate["sd"]
  }
  return (sti.values)
}
