### Unit tests of function likTests

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

test.likTests <- function(){
  ## Poisson Distribution
  data(Polio)
  y <- Polio[, 2]
  X <- as.matrix(Polio[, 3:8])
  glarmamod <- glarma(y, X, thetaLags = c(1,2,5), type = "Poi", method = "FS",
                      residuals = "Pearson", maxit = 100, grad = 2.22e-16)

  lik.tests1 <- likTests(glarmamod)

  ## Test for LR test
  standard.LRTest1 <- c(27.1926023919697, 5.36460455780041e-06)
  checkTrue(abs(lik.tests1[1, 1] - standard.LRTest1[1]) < 10^(-7))
  checkTrue(abs(lik.tests1[1, 2] - standard.LRTest1[2]) < 10^(-7))

  ## Test for Wald test
  standard.WaldTest1 <- c(38.1193252957133, 2.66674532456435e-08)
  checkTrue(abs(lik.tests1[2, 1] - standard.WaldTest1[1]) < 10^(-7))
  checkTrue(abs(lik.tests1[2, 2] - standard.WaldTest1[2]) < 10^(-7))


  ## Binomial Distribution
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
                       na.action = na.omit,x = TRUE)

  X <- glm.LCRobbery$x


  ## Newton Raphson
  glarmamod <- glarma(Y, X, phiLags = c(1), type = "Bin", method = "NR",
                      residuals = "Pearson", maxit = 100, grad = 1e-6)

  lik.tests2 <- likTests(glarmamod)
  standard.LRTest2 <- c(6.11042760345322, 0.0134386612016913)
  checkTrue(abs(lik.tests2[1, 1] - standard.LRTest2[1]) < 10^(-7))
  checkTrue(abs(lik.tests2[1, 2] - standard.LRTest2[2]) < 10^(-7))
  standard.WaldTest2 <- c(6.14430921515927, 0.0131835654128953)
  checkTrue(abs(lik.tests2[2, 1] - standard.WaldTest2[1]) < 10^(-7))
  checkTrue(abs(lik.tests2[2, 2] - standard.WaldTest2[2]) < 10^(-7))

  ## Negative Binomial Distribution
  data(Asthma)
  y <- Asthma[, 1]
  X <- as.matrix(Asthma[, 2:16])
  glarmamod <- glarma(y, X, thetaLags = 7, type = "NegBin", method = "NR",
                      residuals = "Pearson", maxit = 100, grad = 1e-6)

  lik.tests3 <- likTests(glarmamod)
  standard.LRTest3 <- c(7.04717827753939, 0.00793902172903205)
  checkTrue(abs(lik.tests3[1, 1] - standard.LRTest3[1]) < 10^(-7))
  checkTrue(abs(lik.tests3[1, 2] - standard.LRTest3[2]) < 10^(-7))
  standard.WaldTest3 <- c(5.14686152737022, 0.0232884286514632)
  checkTrue(abs(lik.tests2[2, 1] - standard.WaldTest2[1]) < 10^(-7))
  checkTrue(abs(lik.tests2[2, 2] - standard.WaldTest2[2]) < 10^(-7))
}
