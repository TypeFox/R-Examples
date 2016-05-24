# Testing Q-Q and P-P plots
graphicstest.nlFit <- function() {
  param <- c(0, 1, 2, 3)
  names(param) <- c("mu", "sigma", "alpha", "beta")

  # RUnit uses kind = "M-M", normal.kind = "K-R" for RNG. See ?RNGkind
  set.seed(2242, kind = "Marsaglia-Multicarry")
  dataVector <- rnl(1000, param = param) 
  testnlFit <- nlFit(dataVector)
  
  pdf("Q-Q Plot of dataVector.pdf")
  plot.nlFit(testnlFit, which = 3)
  dev.off()

  pdf("P-P Plot of dataVector.pdf")
  plot.nlFit(testnlFit, which = 4)
  dev.off()
}
