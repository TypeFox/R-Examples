### Unit tests of function glarma

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make
test.glarma <- function(){
  ## Purpose: level 1 test of Poisson Distribution with Pearson Residuals.
  data(Polio)
  y <- Polio[, 2]
  X <- as.matrix(Polio[, 3:8])
  glarmamod <- glarma(y, X, thetaLags = c(1, 2, 5), type = "Poi",
                      method = "FS",
                      residuals = "Pearson", maxit = 100, grad = 2.22e-16)
  ARMA.coef <- coef(glarmamod, "ARMA")
  beta.coef <- coef(glarmamod, "beta")
  standard.beta <- c(0.129975397556676, -3.9283713744459, -0.0991261982282867,
                     -0.530844470598217, 0.211127630371824, -0.39323015135847)
  standard.ARMA <- c(0.218459747602161, 0.127231090367536, 0.0872861009061489)
  checkTrue(max(abs(beta.coef - standard.beta)) < 10^(-7))
  checkTrue(max(abs(ARMA.coef - standard.ARMA)) < 10^(-7))

  ## Purpose: level 1 test of Poisson Distribution with Score Residuals
  y <- Polio[, 2]
  X <- as.matrix(Polio[, 3:8])
  glarmamod <- glarma(y, X, thetaLags = c(1, 2, 5), type = "Poi",
                      method = "FS",
                      residuals = "Score", maxit = 100, grad = 2.22e-16)
  ARMA.coef <- coef(glarmamod, "ARMA")
  beta.coef <- coef(glarmamod, "beta")
  standard.beta <- c(0.0437942685983466, -3.89976137445958,
                     -0.007277992577949, -0.588309451810896,
                     0.293551629091528, -0.283751085205294)
  standard.ARMA <- c(0.300327727382111, 0.236693180005699, 0.0182432087875542)
  checkTrue(max(abs(beta.coef - standard.beta)) < 10^(-7))
  checkTrue(max(abs(ARMA.coef - standard.ARMA)) < 10^(-7))

  ## Purpose: level 1 test of Binomial Distribution with Pearson Residuals
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



  ARMA.coef <- coef(glarmamod, "ARMA")
  beta.coef <- coef(glarmamod, "beta")
  standard.beta <- c(-0.274683503145766, 0.822033012059457, -0.35677153900415,
                     -0.50038714171842)
  standard.ARMA <- 0.0817517215109488
  checkTrue(max(abs(beta.coef - standard.beta)) < 10^(-7))
  checkTrue(max(abs(ARMA.coef - standard.ARMA)) < 10^(-7))

  ## Purpose: level 1 test of Negative Binomial Distribution
  ## with Pearson Residuals
  data(Asthma)
  y <- Asthma[, 1]
  X <- as.matrix(Asthma[, 2:16])
  glarmamod <- glarma(y, X, thetaLags = 7, type = "NegBin", method = "NR",
                      residuals = "Pearson", maxit = 100, grad = 1e-6)
  NB.coef <- coef(glarmamod, "NB")
  ARMA.coef <- coef(glarmamod, "ARMA")
  beta.coef <- coef(glarmamod, "beta")
  standard.NB <- 37.1894834266554
  standard.ARMA <- 0.0439191964354301
  standard.beta <- c(0.583971105773007, 0.194554270324514, 0.229989870050125,
                     -0.214500792563189, 0.177283114311301, 0.168433728455995,
                     -0.10403564277627,
                     0.199030077083204, 0.130872744740543, 0.0858677507474468,
                     0.170818285353772,
                     0.252758864156398, 0.305721203236321, 0.436070616091362,
                     0.114120291126706)
  checkTrue(max(abs(NB.coef - standard.NB)) < 10^(-7))
  checkTrue(max(abs(beta.coef - standard.beta)) < 10^(-7))
  checkTrue(max(abs(ARMA.coef - standard.ARMA)) < 10^(-7))


}
