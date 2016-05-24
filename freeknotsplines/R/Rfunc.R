library(splines)

setClass("freekt", 
    representation(x = "numeric", y = "numeric", degree = "integer",
    seed = "integer", stream = "integer", lambda = "numeric",
    optknot = "numeric", tracehat = "numeric", GCV = "numeric",
    GSJS = "numeric", call = "call"))

freelsgen <- function(x, y, degree, numknot, seed = 5, stream = 0)
{
  n <- length(x)
  ord <- degree + 1
  optknot <- rep(0, times = numknot)
  tracehat <- 0
  GCV <- 0
  GSJS <- 0
  result <- .C("freelsgen", as.integer(n), as.double(x), as.double(y), 
               as.integer(ord), as.integer(numknot), as.integer(seed), 
               as.integer(stream), as.double(optknot), 
               as.double(tracehat), as.double(GCV), as.double(GSJS))
  answer <- new("freekt", x = x, y = y, degree = as.integer(degree),
                seed = as.integer(seed), stream = as.integer(stream),
                lambda = 0, optknot = result[[8]], 
                tracehat = result[[9]], GCV = result[[10]], 
                GSJS = result[[11]], call = match.call())
  return(answer)
}

freelsgold <- function(x, y, degree, numknot, seed = 5, stream = 0)
{
  n <- length(x)
  ord <- degree + 1
  optknot <- rep(0, times = numknot)
  tracehat <- 0
  GCV <- 0
  GSJS <- 0
  result <- .C("freelsgold", as.integer(n), as.double(x), as.double(y), 
               as.integer(ord), as.integer(numknot), as.integer(seed), 
               as.integer(stream), as.double(optknot), 
               as.double(tracehat), as.double(GCV), as.double(GSJS))
  answer <- new("freekt", x = x, y = y, degree = as.integer(degree),
                seed = as.integer(seed), stream = as.integer(stream),
                lambda = 0, optknot = result[[8]], 
                tracehat = result[[9]], GCV = result[[10]], 
                GSJS = result[[11]], call = match.call())
  return(answer)
}

freepsgen <- function(x, y, degree, numknot, seed = 5, stream = 0)
{
  n <- length(x)
  ord <- degree + 1
  optknot <- rep(0, times = numknot)
  lambda <- 0
  tracehat <- 0
  GCV <- 0
  GSJS <- 0
  result <- .C("freepsgen", as.integer(n), as.double(x), as.double(y), 
               as.integer(ord), as.integer(numknot), as.integer(seed), 
               as.integer(stream), as.double(lambda), 
               as.double(optknot), as.double(tracehat), 
               as.double(GCV), as.double(GSJS))
  answer <- new("freekt", x = x, y = y, degree = as.integer(degree),
               seed = as.integer(seed), stream = as.integer(stream),
               lambda = result[[8]], optknot = result[[9]], 
               tracehat = result[[10]], GCV = result[[11]], 
               GSJS = result[[12]], call = match.call())
  return(answer)
}

freepsgold <- function(x, y, degree, numknot, seed = 5, stream = 0)
{
  n <- length(x)
  ord <- degree + 1
  optknot <- rep(0, times = numknot)
  lambda <- 0
  tracehat <- 0
  GCV <- 0
  GSJS <- 0
  result <- .C("freepsgold", as.integer(n), as.double(x), as.double(y), 
               as.integer(ord), as.integer(numknot), as.integer(seed), 
               as.integer(stream), as.double(lambda), 
               as.double(optknot), as.double(tracehat), 
               as.double(GCV), as.double(GSJS))
  answer <- new("freekt", x = x, y = y,  degree = as.integer(degree),
               seed = as.integer(seed), stream = as.integer(stream),
               lambda = result[[8]], optknot = result[[9]], 
               tracehat = result[[10]], GCV = result[[11]], 
               GSJS = result[[12]], call = match.call())
  return(answer)
}

chgbasismat <- function(knot, ord)
{
  dimmat <- length(knot) - ord
  answer <- matrix(0, nrow = dimmat, ncol = dimmat)
  for (j in 0:(ord-1))
  {
      brow <- splineDesign(knot, knot[1], ord, j)
      brow <- as.vector(brow/factorial(j))
      answer[j + 1, ] <- brow
  }
  nknot <- sort(-1*knot)
  for (j in 1:(dimmat - ord))
  {
      brow <- splineDesign(knot, knot[ord + j], ord, ord - 1)
      brow2 <- splineDesign(nknot, nknot[length(knot) - ord - (j - 1)],
               ord, ord - 1)
      brow2 <- brow2[dimmat:1]
      brow <- brow + (-1)^ord * brow2
      brow <- as.vector(brow/factorial(ord - 1))
      answer[ord + j, ] <- brow
  }
  return(answer)
}

coef.freekt <- function(object, ...)
{
  xdat <- object@x
  ydat <- object@y
  optknot <- object@optknot
  ord <- object@degree + 1
  lambda <- object@lambda
  fulloptknot <- c(rep(min(xdat), ord), optknot, rep(max(xdat), ord))  # includes endpoints
  Xmat <- splineDesign(fulloptknot, xdat, ord)
  if ((lambda == 0) | (length(optknot) == 0))
    coef <- solve(t(Xmat)%*%Xmat, t(Xmat)%*%ydat) 
  else
    {
      numknots <- length(optknot)
      Amat <- chgbasismat(fulloptknot, ord)
      Istar <- diag(c(rep(0, times = ord), rep(1, times = numknots))) 
      coef <- solve(t(Xmat)%*%Xmat + lambda*t(Amat)%*%Istar%*%Amat, 
                    t(Xmat)%*%ydat) 
    }   
  return(coef)
}

fitted.freekt <- function(object, xfit = object@x, ...)
{        
  xdat <- object@x
  ydat <- object@y
  optknot <- object@optknot
  ord <- object@degree + 1
  fulloptknot <- c(rep(min(xdat), ord), optknot, rep(max(xdat), ord))  # includes endpoints
  coef <- coef.freekt(object)
  yfit <- splineDesign(fulloptknot, xfit, ord) %*%coef
  return(yfit)
}

residuals.freekt <- function(object, ...)
{
  fit <- fitted.freekt(object)
  return(object@y - fit)
}

plot.freekt <- function(x, xfit = x@x, linecolor="blue", lwd = 1, lty = 1, ...)
{
  xfit <- as.vector(xfit)
  yfit <- fitted.freekt(x, xfit)
  plot(x@x, x@y, ...)
  lines(xfit[order(xfit)], yfit[order(xfit)], col=linecolor, 
        lwd = lwd, lty = lty)
}

summary.freekt <- function(object, ...)
{
  answer <- NULL
  if (object@lambda != 0)
  {
     currline <- c(object@lambda, 
        rep(NA, times = length(object@optknot)-1))
     answer <- rbind(answer, currline)
  }
  currline <- object@optknot
  answer <- rbind(answer, currline)
  currline <- c(object@GCV, 
        rep(NA, times = length(object@optknot)-1))
  answer <- rbind(answer, currline)
  RSS <- sum((residuals(object))^2)
  currline <- c(RSS, rep(NA, times = length(object@optknot)-1))
  answer <- rbind(answer, currline)
  if (object@lambda != 0)
       rownames(answer) <- 
           c("Optimal lambda", "Optimal knots", "GCV", "RSS") 
  else
       rownames(answer) <- c("Optimal knots", "GCV", "RSS") 
  colnames(answer) <- rep("", times = length(object@optknot))
  class(answer) <- "table"
  print(answer)
}

AIC.freekt <- function(object, ..., k = 2)
{
  answer <- 0
  RSS <- sum((residuals(object))^2)
  n <- length(object@x)
  npar <- object@tracehat
  answer <- n * log(RSS/n) + k * npar
  return(answer)
}

AICc.freekt <- function(object)
{
  answer <- 0
  RSS <- sum((residuals(object))^2)
  n <- length(object@x)
  npar <- object@tracehat
  answer <- n * log(RSS/n) + 2 * npar + 
            2 * npar * (npar + 1) /(n - npar - 1)
  return(answer)
}

BIC.freekt <- function(object, ...)
{
  answer <- 0
  RSS <- sum((residuals(object))^2)
  n <- length(object@x)
  npar <- object@tracehat
  answer <- n * log(RSS/n) + log(n) * npar
  return(answer)
}

adjGCV.freekt <- function(object, d = 3)
{
  RSS <- sum((residuals(object))^2)
  n <- length(object@x)
  adjtrace <- object@tracehat + d * length(object@optknot)
  answer <- (RSS/n) / (1 - adjtrace / n)^2
  return(answer)     
}

adjAIC.freekt <- function(object)
{
  answer <- 0
  RSS <- sum((residuals(object))^2)
  n <- length(object@x)
  npar <- object@tracehat
  effdim <- 2 * npar - object@degree - 1
  answer <- n * log(RSS/n) +  2 * effdim
  return(answer)
}

fit.search.numknots <- function(x, y, degree, minknot = 1, maxknot = 5, 
                                alg = "LS", search = "genetic",
                                knotnumcrit = "adjGCV", k = 2,
                                d = 3, seed = 5, stream = 0)
{
  bestcrit <- Inf
  funcname <- ""
  answer <- NULL
  if ((alg == "LS") && (search == "genetic"))
    funcname <- "freelsgen"
  if ((alg == "LS") && (search == "golden"))
    funcname <- "freelsgold"
  if ((alg == "PS") && (search == "genetic"))
    funcname <- "freepsgen"
  if ((alg == "PS") && (search == "golden"))
    funcname <- "freepsgold"
  for (numknot in seq(from = minknot, to = maxknot))
  {
      currcall <- call(funcname, x, y, degree, numknot, seed, stream)
      currfit <- eval(currcall)
      currcrit <- switch(knotnumcrit, GCV = currfit$GCV, 
          AIC = AIC(currfit, k = k), adjAIC = adjAIC.freekt(currfit),
          AICc = AICc.freekt(currfit), BIC = BIC(currfit), 
          adjGCV = adjGCV.freekt(currfit, d))
      print(paste("Number of knots = ", numknot, ", ", knotnumcrit, 
          " = ", currcrit, sep = ""), quote = FALSE)
      if (currcrit < bestcrit)
      {
          bestcrit <- currcrit
          answer <- currfit
      }     
  }
  return(answer)
}
