### Graphical Test for Generalized Linear Autocorrelation Moving Average

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make
graphictests.plot.glarma <- function(){
  ## Purpose: Test the plot function with Poisson Distribution
  data(Polio)
  y <- Polio[, 2]
  X <- as.matrix(Polio[, 3:8])
  glarmamod <- glarma(y, X, thetaLags = c(1, 2, 5), type = "Poi", method = "FS",
                      residuals = "Pearson", maxit = 100, grad = 1e-6)
  ## open file for graphical output
  graphicsOutput <- paste(pathReport, "glarma.pdf", sep = "")
  cat("Graphics output in file ", graphicsOutput, "\n")
  pdf(file = graphicsOutput, height = 7, width = 10)
  par(mfrow = c(1, 2))
  plot(glarmamod)

  ## Purpose: Test the plot function with Binomial Distribution
  data(RobberyConvict)
  datalen <- dim(RobberyConvict)[1]
  monthmat <- matrix(0, nrow = datalen, ncol = 12)
  dimnames(monthmat) <- list(NULL, c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  months <- unique(months(strptime(RobberyConvict$Date, "%m/%d/%Y"),
                          abbreviate = TRUE))

  for (j in 1:12) {
    monthmat[months(strptime(RobberyConvict$Date, "%m/%d/%Y"),
                    abbreviate = TRUE) == months[j], j] <-1
  }

  RobberyConvict <- cbind(rep(1, datalen), RobberyConvict, monthmat)
  rm(monthmat)

  ## LOWER COURT ROBBERY
  y1 <- RobberyConvict$LC.Y
  n1 <- RobberyConvict$LC.N

  Y <- cbind(y1, n1 - y1)

  glm.LCRobbery <- glm(Y ~ -1 + Incpt + Step.2001 +
                       I(Feb + Mar + Apr + May + Jun + Jul) +
                       I(Aug + Sep + Oct + Nov + Dec),
                       data = RobberyConvict, family = binomial(link = logit),
                       na.action = na.omit, x = TRUE)

  X <- glm.LCRobbery$x


  ## Newton Raphson
  glarmamod <- glarma(Y, X, phiLags = c(1), type = "Bin", method = "NR",
                      residuals = "Pearson", maxit = 100, grad = 1e-6)
  plot(glarmamod)

  ## Purpose: Test the plot function with Negative Binomial Distribution
  data(Asthma)
  y <- Asthma[, 1]
  X <- as.matrix(Asthma[, 2:16])
  glarmamod <- glarma(y, X, thetaLags = 7, type = "NegBin", method = "NR",
                      residuals = "Pearson", maxit = 100, grad = 1e-6)
  plot(glarmamod)

  dev.off()

  return()
}
