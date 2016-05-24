# Testing summary.nlFit
test.summary.nlFit <- function() {
  param <- c(0, 1, 2, 3)
  names(param) <- c("mu", "sigma", "alpha", "beta")

  # RUnit uses kind = "M-M", normal.kind = "K-R" for RNG. See ?RNGkind
  set.seed(2242, kind = "Marsaglia-Multicarry")
  dataVector <- rnl(1000, param = param)

  testnlFit <- nlFit(dataVector, hessian = TRUE)
  testSummary <- summary(testnlFit)
  checkTrue(is.numeric(testSummary$sds))
}
