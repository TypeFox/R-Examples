# Testing nlFitStart
test.nlFitStart <- function() {
  param <- c(0, 1, 2, 3)
  names(param) <- c("mu", "sigma", "alpha", "beta")

  # RUnit uses kind = "M-M", normal.kind = "K-R" for RNG. See ?RNGkind
  set.seed(2242, kind = "Marsaglia-Multicarry")
  dataVector <- rnl(1000, param = param)

  # Testing user-supplied parameters
  checkEquals(nlFitStart(dataVector, paramStart = param,
                         startValues = "US")$paramStart, param)

  # Testing generated starting parameters, also tests nlFitStartMoM
  testParamStart <- nlFitStart(dataVector)$paramStart
  values <- c(mu = 0.2985418, sigma = 1.0229485,
              alpha = 5.2171185, beta = 3.2115507)
  checkEquals(testParamStart, values)
}
