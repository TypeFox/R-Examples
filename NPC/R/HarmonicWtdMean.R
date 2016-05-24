HarmonicWtdMean <- function (y, tr, tl, block, ...) {
  ## Block-Specific Mean Differences Weighted by Harmonic Mean Sample Size
  df <- data.frame(y, tr, block)
  gdf <- dplyr::group_by(df, block)
  sdf <- dplyr::summarise(gdf,
                          diff = mean.default(y[tr == tl]) -
                              mean.default(y[tr != tl]),
                          Nt <- sum(tr == tl),
                          Nc <- sum(tr != tl),
                          wt = 2 * Nt * Nc / (Nt + Nc))
  return(stats::weighted.mean(sdf$diff, sdf$wt))
}
