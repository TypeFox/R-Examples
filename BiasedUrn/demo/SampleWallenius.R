# SampleWallenius.R
# This demo makes random samples from Wallenius' noncentral hypergeometric 
# distribution and compares measured and expected frequencies

require(BiasedUrn)
require(stats)

MakeSamples <- function(m1, m2, n, odds) {
  nsamp <- 100000                  # Desired number of samples from distribution
  xmin  <- minHypergeo(m1, m2, n)  # Lower limit for x
  xmax  <- maxHypergeo(m1, m2, n)  # Upper limit for x
  
  # Make nsamp samples from Wallenius' distribution
  X <- rWNCHypergeo(nsamp, m1, m2, n, odds)
  
  # Get table of frequencies
  XTab <- as.data.frame(table(X))
  
  # Relative frequencies
  XTab$Freq <- XTab$Freq / nsamp
  
  # Get expected frequencies
  XTab$Expected <- dWNCHypergeo(as.integer(levels(XTab$X)), m1, m2, n, odds)
  
  print("X frequencies in Wallenius' noncentral hypergeometric distribution")
  
  # List measured vs. expected frequencies
  # (How do I get rid of the row names?)
  print(XTab, digits=5)
  
  # Draw histogram
  # (Why does my histogram show densities bigger than 1?)
  hist(X, freq=FALSE)
}

MakeSamples(6, 8, 5, 1.5)