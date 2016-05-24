GAMart <- function(n = 500, sd = 0.1, seed = TRUE)
{
  if(seed) set.seed(111)

  ## help scale function
  hs <- function(x, min = 0.1, max = 0.6) {
    x <- if(length(unique(x)) > 1) {
      (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE)) * (max - min) + min
    } else x
    x
  }
     
  ## (1) regressors
  ## spatial
  n2 <- ceiling(sqrt(n) / 4)
  d <- expand.grid("long" = seq(0, 1, length = n2), "lat" = seq(0, 1, length = n2))
  d$id <- factor(1:nrow(d))
  d <- rep(d, ceiling(n / (n2 * n2)))
  d <- data.frame(
    "long" = unlist(d[grepl("long", names(d))]),
    "lat" = unlist(d[grepl("lat", names(d))]),
    "id" = unlist(d[grepl("id", names(d))])
  )
  d <- d[1:n, ]
  
  ## other covariates
  d$x1 <- runif(n, 0, 1)
  d$x2 <- runif(n, 0, 1)
  d$x3 <- runif(n, 0, 1)
  d$fac <- factor(sample(1:3, size = n, replace = TRUE), labels = c("low", "medium", "high"))

  ## (4) functions
  ## linear
  f1 <- function(x) {
    hs(-1.2 * x, -1, 1)
  }

  ## doublemode
  f2 <- function(x) {
    hs(1.3 * (120 * x * exp(-10 * x) + 2.75 * x^2), -1, 1)
  }

  ## quadratic 
  f3 <- function(x) {
    hs(3.5 * (x - mean(x))^2, -1, 1)
  }

  ## spatial
  f4 <- function(long, lat) {
    hs(sin(hs(long, -3, 3)) * cos(hs(lat, -3, 3)), -1, 1)
  }

  ## random
  f5 <- function(id) {
    hs(rnorm(length(unique(id)), sd = 0.2)[id], -0.2, 0.2)
  }

  ## factor
  f6 <- function(fac) {
    hs(sort(rnorm(length(unique(fac)), sd = 1))[fac], -0.5, 0.5)
  }
  
  ## response
  d$eta <- with(d, hs(f1(x1) + f2(x2) + f3(x3) + f4(long, lat) + f5(id) + f6(fac), -1, 1))
  d$err <- rnorm(n, sd = sd)
  d$num <- with(d, eta + err)
  d$bin <- cut(d$num, quantile(d$num, probs = c(0, 0.5, 1)), labels = c("no", "yes"), include.lowest = TRUE)
  d$cat <- cut(d$num, quantile(d$num), labels = c("none", "low", "medium", "high"), include.lowest = TRUE)
  d <- d[, c("num", "bin", "cat", "eta", "x1", "x2", "x3", "fac", "id", "long", "lat", "err")]

  d
}

