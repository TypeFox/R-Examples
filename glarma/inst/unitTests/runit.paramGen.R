### Unit tests of initial parameters generating functions

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

test.paramGen <- function(){
  ## Poisson Distribution
  data(Polio)
  y <- Polio[, 2]
  X <- as.matrix(Polio[, 3:8])

  ## generate the theta vector
  theta.lags <- c(1, 2, 5)
  theta.init <- c(0, 0, 0)
  theta <- thetaGen(theta.lags, theta.init)

  standard.theta.lags <- c(1, 2, 5)
  standard.theta.init <- structure(c(0, 0, 0),
                                   .Names = c("theta_1", "theta_2", "theta_5"))

  checkEquals(theta[[1]], standard.theta.lags, tolerance = 10^(-7))
  checkEquals(theta[[2]], standard.theta.init, tolerance = 10^(-7))

  ## generate the vector of phi
  phi.lags <- rep(0, 0)
  phi.init <- rep(0, 0)
  phi <- phiGen(phi.lags, phi.init)

  standard.phi.lags <- numeric(0)
  standard.phi.init <- structure(numeric(0), .Names = character(0))

  checkEquals(phi[[1]], standard.phi.lags, tolerance = 10^(-7))
  checkEquals(phi[[2]], standard.phi.init, tolerance = 10^(-7))

  ## generate the delta vector
  delta <- deltaGen(y = y, X = X, phiInit = phi[[2]],
                    thetaInit = theta[[2]], type = "Poi",
                    alpha = 1)

  standard.delta <- structure(c(0.206938270423236, -4.79866147644597,
                                -0.148733251930304, -0.531876816689963,
                                0.169099793073371, -0.432143521504692, 0,
                                0, 0),
                              .Names = c("Intcpt", "Trend", "CosAnnual",
                                         "SinAnnual", "CosSemiAnnual",
                                         "SinSemiAnnual", "theta_1",
                                         "theta_2", "theta_5"))

  checkEquals(delta, standard.delta, tolerance = 10^(-7))
}
