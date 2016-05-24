# Author: P. Poncet

.kernelList <- c("biweight",
                 "chernoff",
                 "cosine",
                 "eddy",
                 "epanechnikov",
                 "gaussian",
                 "optcosine",
                 "rectangular",
                 "triangular",
                 "uniform")

.kernel.biweight <-
function(x,
         ...)
{
  a <- sqrt(7)
  ax <- abs(x)
  k <- ifelse(ax < a, (15/16) * (1 - (ax/a)^2)^2/a, 0)
  structure(list(x = x, k = k, name = "biweight",
            properties = list(integral.K = 1, integral.K2 = 1/2, continuous = Inf, differentiable = Inf)),
            class = ".kernel")
}

# Derivative
.kernel.dbiweight <-
function(x,
         ...)
{
  a <- sqrt(7)
  ax <- abs(x)
  k <- ifelse(ax < a, -(15/4) * x * (1 - (ax/a)^2)/a^3, 0)
  structure(list(x = x, k = k, name = "dbiweight",
            properties = list(integral.K = NA, integral.K2 = 15/(49*a), continuous = Inf, differentiable = Inf)),
            class = ".kernel")
}

.kernel.uniform <- .kernel.chernoff <-
function(x,
         ...)
{
  k <- ifelse(abs(x) <= 1, 1/2, 0)
  structure(list(x = x, k = k, name = "chernoff",
            properties = list(integral.K = 1, integral.K2 = 1/2, continuous = 0, differentiable = 0)),
            class = ".kernel")
}

# Derivative
.kernel.duniform <- .kernel.dchernoff <-
function(x,
         ...)
{
  k <- 0
  structure(list(x = x, k = k, name = "dchernoff",
            properties = list(integral.K = 0, integral.K2 = 0, continuous = Inf, differentiable = Inf)),
            class = ".kernel")
}

.kernel.cosine <-
function(x,
         ...)
{
  a <- 1/sqrt(1/3 - 2/pi^2)
  k <- ifelse(abs(x) < a, (1 + cos(pi*x/a))/(2*a), 0)
  structure(list(x = x, k = k, name = "cosine",
            properties = list(integral.K = 1, integral.K2 = 3/(4*a), continuous = Inf, differentiable = Inf)),
            class = ".kernel")
}

# Derivative                          
.kernel.dcosine <-
function(x,
         ...)
{
  a <- 1/sqrt(1/3 - 2/pi^2)
  k <- ifelse(abs(x) < a, -(pi/(2*a^2))*sin(pi*x/a), 0)
  structure(list(x = x, k = k, name = "dcosine",
            properties = list(integral.K = 0, integral.K2 = pi^2/(4*a^3), continuous = Inf, differentiable = Inf)),
            class = ".kernel")
}

.kernel.eddy <-
function(x,
         ...)
{
  #ax <- abs(x)
  k <- ifelse(abs(x) <= 1, (15/32) * (3 - 10*x^2 + 7*x^4), 0)
  structure(list(x = x, k = k, name = "eddy",
            properties = list(integral.K = 1, integral.K2 = 1.25, continuous = 1, differentiable = 0)),
            class = ".kernel")
}

# Derivative
.kernel.deddy <-
function(x,
         ...)
{
  #ax <- abs(x)
  k <- ifelse(abs(x) <= 1, (15/32) * (-20*x + 28*x^3), 0)
  structure(list(x = x, k = k, name = "deddy",
            properties = list(integral.K = 0, integral.K2 = 9.375, continuous = 0, differentiable = 0)),
            class = ".kernel")
}

.kernel.epanechnikov <-
function(x,
         ...)
{
  #ax <- abs(x)
  k <- ifelse(abs(x) <= 1, (3/4) * (1 - x^2), 0)
  structure(list(x = x, k = k, name = "epanechnikov",
            properties = list(integral.K = 1, integral.K2 = 3/5, continuous = 1, differentiable = 0)),
            class = ".kernel")
}

# Derivative
.kernel.depanechnikov <-
function(x,
         ...)
{
  #ax <- abs(x)
  k <- ifelse(abs(x) <= 1, (-3*x/2), 0)
  structure(list(x = x, k = k, name = "depanechnikov",
            properties = list(integral.K = 0, integral.K2 = 3/2, continuous = 0, differentiable = 0)),
            class = ".kernel")
}

.kernel.gaussian <-
function(x,
         ...)
{
  k <- dnorm(x)
  structure(list(x = x, k = k, name = "gaussian",
            properties = list(integral.K = 1, integral.K2 = 1/(2*sqrt(pi)), continuous = Inf, differentiable = Inf)),
            class = ".kernel")
}

# Derivative
.kernel.dgaussian <-
function(x,
         ...)
{
  k <- -x*dnorm(x)
  structure(list(x = x, k = k, name = "dgaussian",
            properties = list(integral.K = 0, integral.K2 = 0.1410474, continuous = Inf, differentiable = Inf)),
            class = ".kernel")
}

.kernel.optcosine <-
function(x,
         ...)
{
  a <- 1/sqrt(1 - 8/pi^2)
  k <- ifelse(abs(x) < a, (pi/4) * cos(pi * x/(2*a))/a, 0)
  structure(list(x = x, k = k, name = "optcosine",
            properties = list(integral.K = 1, integral.K2 = pi^2/(16*a), continuous = 1, differentiable = 0)),
            class = ".kernel")
}

# Derivative
.kernel.doptcosine <-
function(x,
         ...)
{
  a <- 1/sqrt(1 - 8/pi^2)
  k <- ifelse(abs(x) < a, -(pi^2/(8*a^2)) * sin(pi * x/(2*a)), 0)
  structure(list(x = x, k = k, name = "doptcosine",
            properties = list(integral.K = 0, integral.K2 = pi^4/(64*a^3), continuous = 0, differentiable = 0)),
            class = ".kernel")
}

.kernel.rectangular <-
function(x,
         ...)
{
  a <- sqrt(3)
  k <- ifelse(abs(x) < a, 0.5/a, 0)
  structure(list(x = x, k = k, name = "rectangular",
            properties = list(integral.K = 1, integral.K2 = 1/2, continuous = 0, differentiable = 0)),
            class = ".kernel")
}

# Derivative
.kernel.drectangular <-
function(x,
         ...)
{
  k <- 0
  structure(list(x = x, k = k, name = "drectangular",
            properties = list(integral.K = 0, integral.K2 = 0, continuous = Inf, differentiable = Inf)),
            class = ".kernel")
}

.kernel.triangular <-
function(x,
         ...)
{
  ax <- abs(x)
  k <- ifelse(ax <= 1, (1 - ax), 0)
  structure(list(x = x, k = k, name = "triangular",
            properties = list(integral.K = 1, integral.K2 = 2/3, continuous = 1, differentiable = 0)),
            class = ".kernel")
}

# Derivative
.kernel.dtriangular <-
function(x,
         ...)
{
  ax <- abs(x)
  k <- ifelse(ax <= 1, -sign(x), 0)
  #k <- -sign(x)
  structure(list(x = x, k = k, name = "dtriangular",
            properties = list(integral.K = 0, integral.K2 = 2, continuous = 0, differentiable = 0)),
            class = ".kernel")
}
