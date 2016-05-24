simulated.wikipedia <- function (N=1429649, length=c(100,1000), seed.rng=42) {
  # distributions based on articles of the English Wikipedia (from Wackypedia corpus)
  # extreme articles were deleted (very short, very long, extremely long words, etc.), leaving 1,429,649 items

  # save current RNG state and set new seed (for reproducible data set)
  runif(1) # make sure there is a .Random.seed
  save.seed <- .Random.seed
  set.seed(seed.rng, kind="default")

  # article lengths (number of tokens) with lognormal distribution, 95% in range specified by length argument
  log.mean <- mean(log(length))
  log.sd <- diff(log(length)) / 4
  tokens.log <- rnorm(N, log.mean, log.sd)
  
  # type and token counts have strong linear correlation on logarithmic scale
  types.log <- 0.42 + 0.81 * tokens.log + rnorm(N, 0, 0.1)

  # average word length has very weak, but significant correlation with log token and type counts
  # (for the simulated data, we exaggerate the correlation a little bit)
  avglen <- 4.14 - 0.53 * tokens.log + 0.65 * types.log + rnorm(N, 0, 0.338)

  # round token and type counts, and calculate TTR = token-type-ratio
  tokens <- round(exp(tokens.log))
  types <- round(exp(types.log))
  ttr <- tokens / types

  # combine into final data frame
  WackypediaStats <- data.frame(tokens=tokens, types=types, ttr=ttr, avglen=avglen)

  # restore RNG state, then return generated data set
  .Random.seed <<- save.seed
  WackypediaStats
}
