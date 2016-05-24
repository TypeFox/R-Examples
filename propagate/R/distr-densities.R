############### distribution fitting ###################
## skew-normal distribution, taken from package 'VGAM'
dsn <- function (x, location = 0, scale = 1, shape = 0, log = FALSE) 
{
  zedd <- (x - location)/scale
  loglik <- log(2) + dnorm(zedd, log = TRUE) + pnorm(shape * zedd, log.p = TRUE) - log(scale)
  if (log) loglik else exp(loglik)  
}

## generalized normal distribution,
## taken from PDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
dgnorm <- function(x, alpha = 1, xi = 1, kappa = -0.1) 
{
  1/(exp(log(1 - (kappa * (x - xi))/alpha)^2/(2 * kappa^2)) * (sqrt(2 * pi) * (alpha - x * kappa + kappa * xi))) 
}

## scaled and shifted t-distribution,
dst <- function (x, mean = 0, sd = 1, df = 2) 
{
  dt((x - mean)/sd, df = df)/sd
}

## Gumbel distribution,
## taken from PDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
dgumbel <- function(x, location = 0, scale = 1) {
  z <- (x - location)/scale
  (1/scale) * exp(-z - exp(-z))
}

## Johnson SU distribution,
## taken from PDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
dJSU <- function (x, xi = 0, lambda = 1, gamma = 1, delta = 1) 
{
  z <- (x - xi)/lambda
  delta/(lambda * sqrt(2 * pi) * sqrt(z^2 + 1)) * exp(-0.5 * (gamma + delta * log(z + sqrt(z^2 + 1)))^2)
}

## Johnson SB distribution,
## taken from PDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
dJSB <- function (x, xi = 0, lambda = 1, gamma = 1, delta = 1) 
{
  z <- (x - xi)/lambda
  delta/(lambda * sqrt(2 * pi) * z * (1 - z)) * exp(-0.5 * (gamma + delta * log(z/(1 - z)))^2)
}

## three-parameter weibull distribution with location gamma,
## taken from PDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
dweibull2 <- function(x, location = 0, shape = 1, scale = 1) {
  (shape/scale) * ((x - location)/scale)^(shape - 1) * exp(-((x - location)/scale)^shape)  
}

## four-parameter beta distribution with boundary parameters,
## taken from PDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
dbeta2 <- function(x, alpha1 = 1, alpha2 = 1, a = 0, b = 0) {
  (1/beta(alpha1, alpha2)) * ((x - a)^(alpha1 - 1) * (b - x)^(alpha2 - 1))/(b - a)^(alpha1 + alpha2 - 1)
}

## triangular distribution,
## taken from PDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
dtriang <- function(x, a = 0, b = 1, c = 2) {
  y <- numeric(length(x))
  y[x < a] <- 0
  y[a <= x & x <= b] <- 2 * (x[a <= x & x <= b] - a)/((c - a) * (b - a))
  y[b < x & x <= c] <- 2 * (c - x[b < x & x <= c])/((c - a) * (c - b))
  y[c < x] <- 0  
  return(y)
}

## trapezoidal distribution,
## taken from PDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
dtrap <- function(x, a = 0, b = 1, c = 2, d = 3) {
  y <- numeric(length(x))
  u <- 2/(d + c - b - a)
  y[x < a] <- 0
  y[a <= x & x < b] <- u * (x[a <= x & x < b] - a)/(b - a)
  y[b <= x & x < c] <- u
  y[c <= x & x < d] <- u * (d - x[c <= x & x < d])/(d - c)  
  y[d <= x] <- 0
  return(y)
}

## curvilinear trapezoidal distribution,
## taken from PDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
dctrap <- function(x, a = 0, b = 1, d = 0.1) {
  y <- numeric(length(x))
  mp <- (a + b)/2
  sw <- (b - a)/2
  y[x < a - d] <- 0
  y[a - d <= x & x <= a + d] <- log((sw + d)/(mp - x[a - d <= x & x <= a + d]))
  y[a + d < x & x <= b - d] <- log((sw + d)/(sw - d))
  y[b - d <= x & x <= b + d] <- log((sw + d)/(x[b - d <= x & x <= b + d] - mp))
  y[b + d < x] <- 0  
  return(y)
}

## Generalized Trapezoidal distribution,
## taken from PDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
dgtrap <- function (x, min = 0, mode1 = 1, mode2 = 2, max = 3, n1 = 2, n3 = 2, alpha = 1, log = FALSE) 
{
  y <- numeric(length(x))
  nc <- (2 * n1 * n3) / ((2 * alpha * (mode1 - min) * n3) + ((alpha + 1) * (mode2 - mode1) * n1 * n3) + 
                        (2 * (max - mode2) * n1))
  y[min <= x & x < mode1] <- nc * alpha * ((x[min <= x & x < mode1] - min) / (mode1 - min))^(n1 - 1)
  y[mode1 <= x & x < mode2] <-  nc * (((1 - alpha) * ((x[mode1 <= x & x < mode2] - mode1) / (mode2 - mode1))) + alpha)
  y[mode2 <= x & x <= max] <- nc * ((max - x[mode2 <= x & x <= max]) / (max - mode2))^(n3 - 1)
    
  if (log)  y <- log(y)
  if (any(is.nan(y))) {
    warning("NaN in dtrapezoid")
  }
  else if (any(is.na(y))) {
    warning("NA in dtrapezoid")
  }
  return(y)
}

## Laplacian distribution,
## taken from PDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
dlaplace <- function(x, mean = 0, sigma = 1) {
  1/(sqrt(2) * sigma) * exp(-(sqrt(2) * abs(x - mean))/sigma)
}

## Arcsine distribution,
## taken from PDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
darcsin <- function(x, a = 0, b = 1) { 
  y <- numeric(length(x))  
  y[x <= a] <- 0
  y[a < x & x < b] <- 1/(pi * sqrt((b - x[a < x & x < b]) * (x[a < x & x < b] - a)))
  y[b <= x] <- 0
  return(y)
}

## von Mises distribution,
## taken from PDF in Mathematica's "Ultimate Univariate Probability Distribution Explorer"
dmises <- function(x, mu = 0, kappa = 1) {  
  y <- numeric(length(x))  
  y[x < -pi + mu] <- 0
  y[-pi + mu <= x & x <= pi + mu] <- exp(kappa * cos(x[-pi + mu <= x & x <= pi + mu] - mu))/(2 * pi * besselI(kappa, 0))
  y[pi + mu < x] <- 0
  return(y)
}