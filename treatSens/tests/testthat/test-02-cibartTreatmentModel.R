context("treatSens treatment model argument")

generateData <- function() {
  oldSeed <-  if (exists(".Random.seed")) .Random.seed else NULL
  set.seed(1125)
  n <- 50
  p <- 3

  theta <- 0.5
  beta.z <- runif(p, -1, 0.5)
  beta.x <- runif(p, -3, 3)
  tau <- 4
  sigma <- 1
  
  zeta.z <- 0.5
  zeta.y <- 4
  
  x <- matrix(rnorm(n * (p - 1)), n)
  u.0 <- rbinom(n, 1, theta)
  
  z <- rbinom(n, 1, pnorm(beta.z[1] + x %*% beta.z[-1] + zeta.z * u.0))
  y <- beta.x[1] + x %*% beta.x[-1] + zeta.y * u.0 + tau * z + rnorm(n, 0, sigma)

  x.test <- cbind(x, 1 - z)
  
  if (!is.null(oldSeed)) .Random.seed <- oldSeed
  
  list(y = y, z = as.double(z), x = x, x.test = x.test, zeta.z = zeta.z, zeta.y = zeta.y, theta = theta)
}


test_that("cibart works with various treatmentModel specification methods", {
  data <- generateData()
  control <- treatSens:::cibartControl(5L, 20L, n.thin = 3L, n.thread = 1L)
  
  result <- tryCatch(treatSens:::cibart(data$y, data$z, data$x, data$x.test, data$zeta.y, data$zeta.z, data$theta, "ATE",
                                        probit(family = "normal"), control), error = function(e) e)
  expect_is(result, "list")
  
  result <- tryCatch(treatSens:::cibart(data$y, data$z, data$x, data$x.test, data$zeta.y, data$zeta.z, data$theta, "ATE",
                                        probit, control), error = function(e) e)
  expect_is(result, "list")
  
  result <- tryCatch(treatSens:::cibart(data$y, data$z, data$x, data$x.test, data$zeta.y, data$zeta.z, data$theta, "ATE",
                                        "probit(family = 'normal')", control), error = function(e) e)
  expect_is(result, "list")
  
  result <- tryCatch(treatSens:::cibart(data$y, data$z, data$x, data$x.test, data$zeta.y, data$zeta.z, data$theta, "ATE",
                                        "probit", control), error = function(e) e)
  expect_is(result, "list")
  
  trtModel <- treatSens:::probit(family = "normal") 
  result <- tryCatch(treatSens:::cibart(data$y, data$z, data$x, data$x.test, data$zeta.y, data$zeta.z, data$theta, "ATE",
                                        trtModel, control), error = function(e) e)
  expect_is(result, "list")
  
  expect_error(treatSens:::cibart(data$y, data$z, data$x, data$x.test, data$zeta.y, data$zeta.z, data$theta, "ATE",
                                  notAnObject, control))
})
  
rm(generateData)
